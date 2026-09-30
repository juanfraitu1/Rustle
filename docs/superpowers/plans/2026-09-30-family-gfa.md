# family_gfa Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A standalone binary `family_gfa` that turns one families-stage run (`mcl_families`) into a Bandage-viewable
GFA per O1 family (exon blocks as nodes, each copy a path) and a joint GFA per fusion across families.

**Architecture:** One self-contained file, `src/bin/family_gfa.rs`, in five sections: inputs (parsers), joining
(union-find over exonic bases driven by O1's own PAF CIGARs), compaction (columns to exon-aligned nodes, paths, tags,
invariants), writers (GFA, colours, paths, summary) and fusions (relations to joint graphs). It uses only `std`, `clap`
and `anyhow`; it does not call the `rustle` library. Integration tests in `tests/family_gfa.rs` write tiny synthetic
catalogs and run the binary.

**Tech Stack:** Rust 2021 (cargo 1.93), clap 4 (derive), anyhow; tests use `tempfile` (dev-dependency) and, when on
`PATH`, `vg` 1.x for GFA validation.

**Spec:** `docs/superpowers/specs/2026-09-30-family-gfa-design.md`

## Global Constraints

- Work only in the isolated clone `/mnt/linuxdisk/home/juanfraitu/rustle_m2`, branch `machine2/family-gfa`. Never
  touch `/mnt/c/Users/jfris/Desktop/Rustle` (another session is editing it).
- New files only: `src/bin/family_gfa.rs`, `tests/family_gfa.rs`, `docs/FAMILY_GFA_SMOKE_2026-09-30.md`. The one
  modification allowed is appending a `[[bin]]` block to `Cargo.toml` (the repo sets `autobins = false`).
- Never stage `bench/soto/soto_replication.py`, `docs/HANDOFF_SECOND_MACHINE_2026-09-30.md` or
  `docs/SOTO_VS_OURS_MEETING_2026-09-30.md` (another session's uncommitted work, copied into this clone). Stage files
  by explicit path, never `git add -A` / `git add .`.
- Build and test with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2`, always `--release`, through
  `tools/rlock.sh heavy`, with cargo output captured to a file under `/mnt/linuxdisk/tmp/`. **Never** use
  `/mnt/linuxdisk/home/juanfraitu/rustle_target`: that target dir belongs to the other session and building there
  overwrites its binaries.
- Output must be byte-identical across reruns: iterate `BTreeMap`/`BTreeSet` or sorted vectors when writing; use
  `HashMap`/`HashSet` only for lookups.
- GFA 1.0. First two lines of a family GFA, verbatim:
  `H\tVN:Z:1.0` then `# graph view of O1 family <id>, built by family_gfa from O1's own alignments; not used to define families`.
- Wording in docs and messages: "a graph view of an O1 family". Never "a family is one VG".
- Colours (Bandage CSV, header `Node,Colour`): private node (1 copy) `#eb6834`; shared nodes a ramp from `#86b6ef`
  (2 copies) to `#104281` (all copies; a 2-copy family's shared nodes are `#104281`); fusion-link end nodes `#e34948`.
- Commit messages end with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
  ```
  Commit with the main repo's identity:
  `git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit ...`

## Review Focus

1. A PAF from a different run whose coordinates run past a locus span: expect an error naming both loci, not a panic
   (test in Task 2).
2. Soft-masked (lowercase) bases in `.loci.fa`: expect upper-casing on read, so identities are not deflated (test in
   Task 1).
3. `--families` naming a family that is not in `.clusters.tsv`: expect an error naming it (test in Task 1).
4. A relation unit whose locus is absent from the loci files: expect a `reason` in `fusions.tsv`, not a crash
   (test in Task 5).
5. A relation whose `unit_families` disagrees with `.clusters.tsv` (6 of 349 units in the real chr16 run): expect
   `.clusters.tsv` to win and the disagreement written to the `note` column (test in Task 5).

## Common commands

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2
export CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2
# unit tests (in the binary)
tools/rlock.sh heavy cargo test --release --bin family_gfa > /mnt/linuxdisk/tmp/fg_unit.log 2>&1; echo exit=$?; tail -25 /mnt/linuxdisk/tmp/fg_unit.log
# integration tests
tools/rlock.sh heavy cargo test --release --test family_gfa > /mnt/linuxdisk/tmp/fg_int.log 2>&1; echo exit=$?; tail -25 /mnt/linuxdisk/tmp/fg_int.log
```
The first build compiles the whole crate with fat LTO (several minutes); later ones relink only.

---

### Task 1: Inputs, CLI and load

**Files:**
- Create: `src/bin/family_gfa.rs`
- Modify: `Cargo.toml` (append a `[[bin]]` block at the end)
- Create: `tests/family_gfa.rs`

**Interfaces:**
- Produces (used by every later task): `struct Key { chrom: String, start: u64, end: u64 }` with `Key::parse(&str) -> Result<Key>`,
  `Key::name(&self) -> String`, `Key::len(&self) -> usize`; `struct Locus { key, gene_id: String, strand: char, exons: Vec<(u64, u64)>, seq: Vec<u8> }`
  with `exonic() -> Vec<(usize, usize)>`, `exonic_bp() -> usize`, `transcript() -> Vec<u8>`, `path_name() -> String`;
  `fn comp(u8) -> u8`; `struct Paf { q: Key, qs: usize, qe: usize, rev: bool, t: Key, ts: usize, te: usize, cigar: Vec<(u32, u8)> }`;
  `fn parse_cigar(&str) -> Result<Vec<(u32, u8)>>`;
  `struct Inputs { all: BTreeMap<String, Vec<Key>>, selected: Vec<String>, member_of: HashMap<Key, String>, loci: HashMap<Key, Locus>, pafs: Vec<Paf>, no_cigar: usize }`;
  `fn load(prefix: &str, select: &[String], extra: &[Key]) -> Result<Inputs>` (`extra` = relation unit keys, loaded
  when present, skipped silently when absent).
- Test fixture (in `tests/family_gfa.rs`, reused by Tasks 4 and 5): `fixture() -> Vec<Loc>`, `write_inputs(dir, &[Loc], &[String]) -> String`
  (returns the prefix), `all_pafs(&[Loc]) -> Vec<String>`, `run(&[&str]) -> Output`.

- [ ] **Step 1: Append the bin to Cargo.toml**

Append at the end of `Cargo.toml`:

```toml

[[bin]]
name = "family_gfa"
path = "src/bin/family_gfa.rs"
```

- [ ] **Step 2: Write the file with the input section and failing unit tests**

Create `src/bin/family_gfa.rs`:

```rust
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
```

- [ ] **Step 3: Write the integration test fixture and the load tests**

Create `tests/family_gfa.rs`:

```rust
//! `family_gfa` end to end on a synthetic families-stage run: FAM1 = A (+), B (+, one SNP in exon 1), C (-, the
//! reverse complement of A); FAM2 = D, E identical; F (chr6) is in no family.

use std::collections::HashSet;
use std::path::Path;
use std::process::{Command, Output};

/// Deterministic pseudo-random DNA.
fn dna(n: usize, seed: u64) -> String {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize] as char
        })
        .collect()
}

fn rc(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            o => o,
        })
        .collect()
}

pub struct Loc {
    chrom: &'static str,
    start: u64,
    strand: char,
    /// exons as 1-based inclusive offsets into `seq`
    exons: Vec<(u64, u64)>,
    seq: String,
    gid: &'static str,
    /// family id, or "-" for a locus in no family
    fam: &'static str,
}

impl Loc {
    fn key(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.start + self.seq.len() as u64 - 1)
    }
}

fn fixture() -> Vec<Loc> {
    let s = dna(60, 7);
    let mut b = s.clone().into_bytes();
    b[4] = if b[4] == b'A' { b'C' } else { b'A' };
    let b = String::from_utf8(b).unwrap();
    let t = dna(30, 11);
    vec![
        Loc { chrom: "chr1", start: 1001, strand: '+', exons: vec![(1, 20), (41, 60)], seq: s.clone(), gid: "GA", fam: "FAM1" },
        Loc { chrom: "chr2", start: 2001, strand: '+', exons: vec![(1, 20), (41, 60)], seq: b, gid: "GB", fam: "FAM1" },
        Loc { chrom: "chr3", start: 3001, strand: '-', exons: vec![(1, 20), (41, 60)], seq: rc(&s), gid: "GC", fam: "FAM1" },
        Loc { chrom: "chr4", start: 4001, strand: '+', exons: vec![(1, 30)], seq: t.clone(), gid: "GD", fam: "FAM2" },
        Loc { chrom: "chr5", start: 5001, strand: '+', exons: vec![(1, 30)], seq: t, gid: "GE", fam: "FAM2" },
        Loc { chrom: "chr6", start: 6001, strand: '+', exons: vec![(1, 25)], seq: dna(25, 13), gid: "GF", fam: "-" },
    ]
}

/// Writes PREFIX.{clusters.tsv, loci.gff3, loci.fa, loci.paf} into `dir`; returns PREFIX.
fn write_inputs(dir: &Path, locs: &[Loc], paf: &[String]) -> String {
    let p = dir.join("t.fam").to_string_lossy().to_string();
    let mut cl = String::from("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n");
    let mut g = String::from("##gff-version 3\n");
    let mut fa = String::new();
    for l in locs {
        let end = l.start + l.seq.len() as u64 - 1;
        if l.fam != "-" {
            cl += &format!("{}\t0\t1\t0\tNA\t{}\t{}\t{end}\n", l.fam, l.chrom, l.start);
        }
        g += &format!("{}\t.\tgene\t{}\t{end}\t.\t{}\t.\tID=gene-{};Name={}\n", l.chrom, l.start, l.strand, l.gid, l.gid);
        for &(a, b) in &l.exons {
            g += &format!(
                "{}\t.\texon\t{}\t{}\t.\t{}\t.\tParent=gene-{};gene={}\n",
                l.chrom,
                l.start + a - 1,
                l.start + b - 1,
                l.strand,
                l.gid,
                l.gid
            );
        }
        fa += &format!(">{}\n{}\n", l.key(), l.seq.to_lowercase());
    }
    std::fs::write(format!("{p}.clusters.tsv"), cl).unwrap();
    std::fs::write(format!("{p}.loci.gff3"), g).unwrap();
    std::fs::write(format!("{p}.loci.fa"), fa).unwrap();
    std::fs::write(format!("{p}.loci.paf"), paf.join("")).unwrap();
    p
}

/// One full-length PAF line between two loci of equal length.
fn paf_full(q: &Loc, t: &Loc, strand: char) -> String {
    let n = q.seq.len();
    format!("{}\t{n}\t0\t{n}\t{strand}\t{}\t{n}\t0\t{n}\t{n}\t{n}\t60\tNM:i:0\tcg:Z:{n}M\n", q.key(), t.key())
}

/// All-vs-all records within FAM1 (C against A and B is reverse-strand) and FAM2, in both directions.
fn all_pafs(l: &[Loc]) -> Vec<String> {
    vec![
        paf_full(&l[0], &l[1], '+'),
        paf_full(&l[1], &l[0], '+'),
        paf_full(&l[0], &l[2], '-'),
        paf_full(&l[2], &l[0], '-'),
        paf_full(&l[1], &l[2], '-'),
        paf_full(&l[2], &l[1], '-'),
        paf_full(&l[3], &l[4], '+'),
        paf_full(&l[4], &l[3], '+'),
    ]
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_family_gfa")).args(args).output().unwrap()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

#[test]
fn loads_a_valid_run() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let out = d.path().join("out");
    let o = run(&["--prefix", &p, "--outdir", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(stderr(&o).contains("2 families selected, 5 loci, 8 alignment records"), "{}", stderr(&o));
}

#[test]
fn a_member_missing_from_the_gff3_is_an_error_naming_it() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let mut cl = std::fs::read_to_string(format!("{p}.clusters.tsv")).unwrap();
    cl += "FAM3\t0\t1\t0\tNA\tchr9\t1\t10\n";
    std::fs::write(format!("{p}.clusters.tsv"), cl).unwrap();
    let o = run(&["--prefix", &p, "--outdir", d.path().join("o").to_str().unwrap()]);
    assert!(!o.status.success());
    assert!(stderr(&o).contains("chr9:1-10") && stderr(&o).contains("FAM3"), "{}", stderr(&o));
}

#[test]
fn an_unknown_family_is_an_error_naming_it() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let o = run(&["--prefix", &p, "--outdir", d.path().join("o").to_str().unwrap(), "--families", "FAM1,NOPE"]);
    assert!(!o.status.success());
    assert!(stderr(&o).contains("family NOPE is not in"), "{}", stderr(&o));
}

#[allow(dead_code)]
fn unused_helpers_kept_for_later_tasks(_: HashSet<String>) {}
```

Note: `write_inputs` writes the FASTA in lower case on purpose (Review Focus 2); the unit test
`fasta_uppercases_and_keeps_only_wanted` covers upper-casing directly.

- [ ] **Step 4: Run the tests**

Run the unit and integration commands from "Common commands". Expected: all 5 unit tests and 3 integration tests
PASS. If a test fails, fix the code in `family_gfa.rs` (not the test) until it passes.

- [ ] **Step 5: Commit**

```bash
git add Cargo.toml src/bin/family_gfa.rs tests/family_gfa.rs
git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -F - <<'EOF'
family_gfa: inputs, CLI and load (graph view of O1 families, task 1/6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
EOF
```

---

### Task 2: Joining exonic bases into columns

**Files:**
- Modify: `src/bin/family_gfa.rs` (add a `// Joining` section after the inputs section; add unit tests to `mod tests`)

**Interfaces:**
- Consumes: `Locus`, `Paf`, `Key`, `parse_cigar` from Task 1.
- Produces:
  `#[derive(Clone, Debug, Default, PartialEq)] struct Stats { records: usize, loops: usize, inversion_conflicts: usize, same_dna_bp: usize }`;
  `struct Joined { members: Vec<Locus>, exonic: Vec<Vec<(usize, usize)>>, vid: Vec<Vec<u32>>, owner: Vec<(usize, usize)>, col: Vec<usize>, orient: Vec<u8>, rep: Vec<u32>, same_dna: Vec<bool>, stats: Stats }`;
  `fn join(name: &str, members: Vec<Locus>, pafs: &[&Paf]) -> Result<Joined>` (members sorted by key inside);
  `fn columns(r: &Paf) -> Result<Vec<(usize, usize)>>`.

- [ ] **Step 1: Write the failing unit tests**

Add to `mod tests` in `src/bin/family_gfa.rs`:

```rust
    /// A locus on chrom `c` starting at `start`, exons as 1-based offsets into `seq` (inclusive).
    fn locus(c: &str, start: u64, strand: char, exons: &[(u64, u64)], seq: &str, gid: &str) -> Locus {
        let key = Key { chrom: c.into(), start, end: start + seq.len() as u64 - 1 };
        Locus {
            key,
            gene_id: gid.into(),
            strand,
            exons: exons.iter().map(|&(a, b)| (start + a - 1, start + b - 1)).collect(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn paf(q: &Locus, qs: usize, qe: usize, rev: bool, t: &Locus, ts: usize, te: usize, cigar: &str) -> Paf {
        Paf { q: q.key.clone(), qs, qe, rev, t: t.key.clone(), ts, te, cigar: parse_cigar(cigar).unwrap() }
    }

    fn rc(s: &str) -> String {
        s.bytes().rev().map(|b| comp(b) as char).collect()
    }

    const S16: &str = "GATTACACCGTAGGCT";

    #[test]
    fn dsu_tracks_orientation_and_refuses_inversions() {
        let mut d = Dsu::new(3);
        assert!(d.union(0, 1, 1));
        assert!(d.union(1, 2, 1));
        let (r0, p0) = d.find(0);
        let (r2, p2) = d.find(2);
        assert_eq!(r0, r2);
        assert_eq!(p0 ^ p2, 0); // reversed twice = same orientation
        assert!(!d.union(0, 2, 1));
        assert!(d.union(0, 2, 0));
    }

    #[test]
    fn columns_walk_plus_and_minus_records() {
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 16)], S16, "B");
        assert_eq!(columns(&paf(&a, 0, 5, false, &b, 10, 14, "2M1I2M")).unwrap(), vec![(0, 10), (1, 11), (3, 12), (4, 13)]);
        assert_eq!(columns(&paf(&a, 0, 4, true, &b, 0, 4, "4M")).unwrap(), vec![(3, 0), (2, 1), (1, 2), (0, 3)]);
        assert!(columns(&paf(&a, 0, 4, false, &b, 0, 4, "5M")).is_err());
    }

    #[test]
    fn identical_copies_share_every_column() {
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 16)], S16, "B");
        let r = paf(&a, 0, 16, false, &b, 0, 16, "16M");
        let j = join("F", vec![a, b], &[&r]).unwrap();
        assert_eq!(j.rep.len(), 16);
        assert_eq!(j.stats, Stats { records: 1, loops: 0, inversion_conflicts: 0, same_dna_bp: 0 });
    }

    #[test]
    fn an_alignment_only_in_introns_joins_nothing() {
        let a = locus("c1", 1, '+', &[(1, 4), (13, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 4), (13, 16)], S16, "B");
        let r = paf(&a, 4, 12, false, &b, 4, 12, "8M");
        let j = join("F", vec![a, b], &[&r]).unwrap();
        assert_eq!(j.rep.len(), 16); // 8 exonic bases per copy, none joined
    }

    #[test]
    fn an_internal_repeat_is_counted_as_a_loop() {
        // A carries the same 6 bp in both exons; B aligns to both
        let a = locus("c1", 1, '+', &[(1, 6), (11, 16)], "GATTACCCCCGATTAC", "A");
        let b = locus("c2", 1, '+', &[(1, 6)], "GATTAC", "B");
        let r1 = paf(&b, 0, 6, false, &a, 0, 6, "6M");
        let r2 = paf(&b, 0, 6, false, &a, 10, 16, "6M");
        let j = join("F", vec![a, b], &[&r1, &r2]).unwrap();
        assert_eq!(j.rep.len(), 6);
        assert_eq!(j.stats.loops, 6);
    }

    #[test]
    fn same_dna_joins_are_counted() {
        // two overlapping loci on one chromosome, aligned over their 10 bp overlap
        let s = "ACGTACGTACGGATTACAGG";
        let a = locus("c1", 1, '+', &[(1, 20)], s, "A");
        let b = locus("c1", 11, '+', &[(1, 20)], &format!("{}TTTTTTTTTT", &s[10..]), "B");
        let r = paf(&a, 10, 20, false, &b, 0, 10, "10M");
        let j = join("F", vec![a, b], &[&r]).unwrap();
        assert_eq!(j.stats.same_dna_bp, 10);
        assert_eq!(j.same_dna.iter().filter(|&&x| x).count(), 10);
    }

    #[test]
    fn contradictory_strands_are_inversion_conflicts() {
        // a forward and a reverse record over 5 bp: the middle base is joined to itself both ways
        let a = locus("c1", 1, '+', &[(1, 5)], "ACGTA", "A");
        let b = locus("c2", 1, '+', &[(1, 5)], "ACGTA", "B");
        let r1 = paf(&a, 0, 5, false, &b, 0, 5, "5M");
        let r2 = paf(&a, 0, 5, true, &b, 0, 5, "5M");
        let j = join("F", vec![a, b], &[&r1, &r2]).unwrap();
        assert_eq!(j.stats.inversion_conflicts, 1);
    }

    #[test]
    fn a_record_past_a_locus_span_is_an_error_naming_both_loci() {
        let a = locus("c1", 1, '+', &[(1, 4)], "ACGT", "A");
        let b = locus("c2", 1, '+', &[(1, 4)], "ACGT", "B");
        let r = paf(&a, 0, 6, false, &b, 0, 6, "6M");
        let e = join("F", vec![a, b], &[&r]).err().unwrap().to_string();
        assert!(e.contains("c1:1-4") && e.contains("c2:1-4"), "{e}");
    }

    #[test]
    fn reverse_records_join_with_reversed_orientation() {
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '-', &[(1, 16)], &rc(S16), "B");
        let r = paf(&a, 0, 16, true, &b, 0, 16, "16M");
        let j = join("F", vec![a, b], &[&r]).unwrap();
        assert_eq!(j.rep.len(), 16);
        // every B vertex is reversed relative to its column's representative (an A vertex)
        let b_verts: Vec<usize> = (16..32).collect();
        assert!(b_verts.iter().all(|&v| j.orient[v] == 1));
    }
```

- [ ] **Step 2: Run the unit tests and confirm they fail to compile**

Run the unit test command. Expected: compile errors `cannot find type Dsu`, `cannot find function join`,
`cannot find function columns`.

- [ ] **Step 3: Implement the joining section**

Insert after `read_paf`'s section (before `// Command line`) in `src/bin/family_gfa.rs`:

```rust
// ------------------------------------------------------------------------------------------------------------------
// Joining
// ------------------------------------------------------------------------------------------------------------------

/// Union-find over vertices with an orientation parity: `parity[v]` is whether v is reversed relative to its parent.
struct Dsu {
    parent: Vec<u32>,
    parity: Vec<u8>,
    rank: Vec<u8>,
}

impl Dsu {
    fn new(n: usize) -> Dsu {
        Dsu { parent: (0..n as u32).collect(), parity: vec![0; n], rank: vec![0; n] }
    }

    /// (root, parity of v relative to the root), compressing the path.
    fn find(&mut self, v: u32) -> (u32, u8) {
        let mut path = Vec::new();
        let mut x = v;
        while self.parent[x as usize] != x {
            path.push(x);
            x = self.parent[x as usize];
        }
        let root = x;
        let mut acc = 0u8;
        for &y in path.iter().rev() {
            acc ^= self.parity[y as usize];
            self.parity[y as usize] = acc;
            self.parent[y as usize] = root;
        }
        (root, if path.is_empty() { 0 } else { self.parity[v as usize] })
    }

    /// Joins a and b so that orientation(a) = orientation(b) XOR p. Returns false, changing nothing, when they are
    /// already joined with the other parity: an inversion conflict.
    fn union(&mut self, a: u32, b: u32, p: u8) -> bool {
        let (ra, pa) = self.find(a);
        let (rb, pb) = self.find(b);
        if ra == rb {
            return pa ^ pb == p;
        }
        let x = p ^ pa ^ pb;
        let (hi, lo) = if self.rank[ra as usize] >= self.rank[rb as usize] { (ra, rb) } else { (rb, ra) };
        self.parent[lo as usize] = hi;
        self.parity[lo as usize] = x;
        if self.rank[ra as usize] == self.rank[rb as usize] {
            self.rank[hi as usize] += 1;
        }
        true
    }
}

/// Aligned column pairs (query offset, target offset) of one record, offsets 0-based into each locus span. `M`/`=`/`X`
/// pair bases, `I` consumes the query, `D` the target; a `-` record walks the query from its end backwards. The CIGAR
/// must consume exactly the record's two intervals.
fn columns(r: &Paf) -> Result<Vec<(usize, usize)>> {
    let (qlen, tlen) = r.cigar.iter().fold((0usize, 0usize), |(q, t), &(n, op)| match op {
        b'M' | b'=' | b'X' => (q + n as usize, t + n as usize),
        b'I' => (q + n as usize, t),
        b'D' | b'N' => (q, t + n as usize),
        _ => (q, t),
    });
    if r.qe < r.qs || r.te < r.ts || qlen != r.qe - r.qs || tlen != r.te - r.ts {
        bail!(
            "alignment {} vs {}: CIGAR covers {qlen}/{tlen} bp but the record spans {}/{}",
            r.q.name(),
            r.t.name(),
            r.qe.saturating_sub(r.qs),
            r.te.saturating_sub(r.ts)
        );
    }
    let mut out = Vec::new();
    let (mut qo, mut t) = (0usize, r.ts);
    for &(n, op) in &r.cigar {
        let n = n as usize;
        match op {
            b'M' | b'=' | b'X' => {
                for k in 0..n {
                    let q = if r.rev { r.qe - 1 - (qo + k) } else { r.qs + qo + k };
                    out.push((q, t + k));
                }
                qo += n;
                t += n;
            }
            b'I' => qo += n,
            b'D' | b'N' => t += n,
            _ => {}
        }
    }
    Ok(out)
}

/// Per-family diagnostics for summary.tsv.
#[derive(Clone, Debug, Default, PartialEq)]
struct Stats {
    /// PAF records between two members
    records: usize,
    /// columns holding two positions of one copy (an internal repeat)
    loops: usize,
    /// joins refused because they would reverse a base against itself
    inversion_conflicts: usize,
    /// columns whose bases all sit at one chromosome position (overlapping loci: the same DNA, not paralogy)
    same_dna_bp: usize,
}

/// A family's exonic bases joined into columns.
struct Joined {
    /// members sorted by key
    members: Vec<Locus>,
    /// per member: (offset, exon index) of each exonic base, genomic order
    exonic: Vec<Vec<(usize, usize)>>,
    /// per member: vertex id at each span offset (u32::MAX = not exonic)
    vid: Vec<Vec<u32>>,
    /// per vertex: (member, offset)
    owner: Vec<(usize, usize)>,
    /// per vertex: its column
    col: Vec<usize>,
    /// per vertex: 1 if reversed relative to its column's representative
    orient: Vec<u8>,
    /// per column: the representative vertex (lowest member, then lowest offset); it fixes the column's orientation
    rep: Vec<u32>,
    /// per column: all its bases at one chromosome position
    same_dna: Vec<bool>,
    stats: Stats,
}

fn join(name: &str, mut members: Vec<Locus>, pafs: &[&Paf]) -> Result<Joined> {
    members.sort_by(|a, b| a.key.cmp(&b.key));
    let index: HashMap<Key, usize> = members.iter().enumerate().map(|(i, m)| (m.key.clone(), i)).collect();
    let exonic: Vec<Vec<(usize, usize)>> = members.iter().map(|m| m.exonic()).collect();
    let mut vid = Vec::with_capacity(members.len());
    let mut owner: Vec<(usize, usize)> = Vec::new();
    for (i, m) in members.iter().enumerate() {
        let mut v = vec![u32::MAX; m.seq.len()];
        for &(o, _) in &exonic[i] {
            v[o] = owner.len() as u32;
            owner.push((i, o));
        }
        vid.push(v);
    }
    let n = owner.len();
    let mut dsu = Dsu::new(n);
    let mut stats = Stats::default();
    for r in pafs {
        let (Some(&qi), Some(&ti)) = (index.get(&r.q), index.get(&r.t)) else { continue };
        if qi == ti {
            continue;
        }
        stats.records += 1;
        for (qo, to) in columns(r)? {
            let (Some(&a), Some(&b)) = (vid[qi].get(qo), vid[ti].get(to)) else {
                bail!(
                    "family {name}: alignment {} vs {} runs past a locus span (a PAF from another run?)",
                    r.q.name(),
                    r.t.name()
                );
            };
            if a == u32::MAX || b == u32::MAX {
                continue;
            }
            if !dsu.union(a, b, r.rev as u8) {
                stats.inversion_conflicts += 1;
            }
        }
    }
    // columns numbered by their first vertex (member order, then genomic order): that vertex is the representative
    let mut col_of_root: HashMap<u32, usize> = HashMap::new();
    let mut col = vec![0usize; n];
    let mut orient = vec![0u8; n];
    let mut rep: Vec<u32> = Vec::new();
    let mut rep_par: Vec<u8> = Vec::new();
    for v in 0..n as u32 {
        let (r, p) = dsu.find(v);
        let c = *col_of_root.entry(r).or_insert_with(|| {
            rep.push(v);
            rep_par.push(p);
            rep.len() - 1
        });
        col[v as usize] = c;
        orient[v as usize] = p ^ rep_par[c];
    }
    let ncol = rep.len();
    let mut col_verts: Vec<Vec<u32>> = vec![Vec::new(); ncol];
    for v in 0..n {
        col_verts[col[v]].push(v as u32);
    }
    let gpos = |v: u32| {
        let (m, o) = owner[v as usize];
        (members[m].key.chrom.as_str(), members[m].key.start + o as u64)
    };
    let mut same_dna = vec![false; ncol];
    for (c, vs) in col_verts.iter().enumerate() {
        let mut seen = HashSet::new();
        if vs.iter().any(|&v| !seen.insert(owner[v as usize].0)) {
            stats.loops += 1;
        }
        if vs.len() >= 2 && vs.iter().all(|&v| gpos(v) == gpos(vs[0])) {
            same_dna[c] = true;
            stats.same_dna_bp += 1;
        }
    }
    Ok(Joined { members, exonic, vid, owner, col, orient, rep, same_dna, stats })
}
```

- [ ] **Step 4: Run the unit tests**

Expected: the 9 new tests and the 5 from Task 1 PASS (dead-code warnings for `Joined` fields not yet read are fine).
Run the integration tests too: the 3 from Task 1 still PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bin/family_gfa.rs
git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -F - <<'EOF'
family_gfa: join exonic bases into columns from O1's PAF (task 2/6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
EOF
```

---

### Task 3: Compaction into exon-aligned nodes, paths, tags and invariants

**Files:**
- Modify: `src/bin/family_gfa.rs` (add a `// Compaction` section after `// Joining`; add unit tests)

**Interfaces:**
- Consumes: `join`, `Joined`, `Stats`, `comp`, `Locus` from Tasks 1-2.
- Produces:
  `struct Node { seq: Vec<u8>, copies: usize, traversals: usize, id_mean: f64, id_min: f64, same_dna: bool }`;
  `type Link = (usize, bool, usize, bool);` `fn canon(a: usize, ar: bool, b: usize, br: bool) -> Link`;
  `struct Graph { name: String, members: Vec<Locus>, nodes: Vec<Node>, paths: Vec<Vec<(usize, bool)>>, links: BTreeSet<Link>, mism: Vec<usize>, stats: Stats }`
  (node indices are 0-based; node `n` is written as GFA id `n + 1`; `paths[i]` is member `i`'s path as (node, reversed));
  `fn build(name: &str, members: Vec<Locus>, pafs: &[&Paf]) -> Result<Graph>`; `fn check(g: &Graph) -> Result<()>`;
  `fn oriented(seq: &[u8], rev: bool) -> Vec<u8>`.

- [ ] **Step 1: Write the failing unit tests**

Add to `mod tests`:

```rust
    fn graph(ms: Vec<Locus>, ps: Vec<Paf>) -> Graph {
        let refs: Vec<&Paf> = ps.iter().collect();
        let g = build("F", ms, &refs).unwrap();
        check(&g).unwrap();
        g
    }

    #[test]
    fn identical_copies_give_one_shared_node() {
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 16)], S16, "B");
        let r = paf(&a, 0, 16, false, &b, 0, 16, "16M");
        let g = graph(vec![a, b], vec![r]);
        assert_eq!(g.nodes.len(), 1);
        assert_eq!((g.nodes[0].copies, g.nodes[0].id_mean, g.nodes[0].id_min), (2, 1.0, 1.0));
        assert_eq!(g.paths, vec![vec![(0, false)], vec![(0, false)]]);
    }

    #[test]
    fn a_snp_lowers_identity_without_a_bubble() {
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 16)], "GATTCCACCGTAGGCT", "B"); // A->C at offset 4
        let r = paf(&a, 0, 16, false, &b, 0, 16, "16M");
        let g = graph(vec![a, b], vec![r]);
        assert_eq!(g.nodes.len(), 1);
        assert_eq!(g.nodes[0].id_mean, 15.0 / 16.0);
        assert_eq!(g.mism, vec![0, 1]);
    }

    #[test]
    fn nodes_are_cut_at_exon_boundaries() {
        let a = locus("c1", 1, '+', &[(1, 6), (11, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 6), (11, 16)], S16, "B");
        let r = paf(&a, 0, 16, false, &b, 0, 16, "16M");
        let g = graph(vec![a, b], vec![r]);
        assert_eq!(g.nodes.len(), 2);
        assert_eq!(g.paths[0], vec![(0, false), (1, false)]);
    }

    #[test]
    fn an_extra_exon_is_a_private_node() {
        let s = "AAAACCCCGGGGTTTTACGT";
        let a = locus("c1", 1, '+', &[(1, 4), (9, 12), (17, 20)], s, "A");
        let b = locus("c2", 1, '+', &[(1, 4), (9, 12), (13, 16), (17, 20)], s, "B");
        let r = paf(&a, 0, 20, false, &b, 0, 20, "20M");
        let g = graph(vec![a, b], vec![r]);
        assert_eq!(g.nodes.len(), 4);
        let private: Vec<&Node> = g.nodes.iter().filter(|n| n.copies == 1).collect();
        assert_eq!(private.len(), 1);
        assert_eq!(private[0].seq, b"TTTT".to_vec());
        assert_eq!(g.paths[0].len(), 3);
        assert_eq!(g.paths[1].len(), 4);
    }

    #[test]
    fn a_reverse_strand_copy_spells_its_transcript() {
        // B is A's reverse complement on the minus strand: same transcript, walked forward through A's node
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '-', &[(1, 16)], &rc(S16), "B");
        let r = paf(&a, 0, 16, true, &b, 0, 16, "16M");
        let g = graph(vec![a, b], vec![r]);
        assert_eq!(g.nodes.len(), 1);
        assert_eq!(g.paths[1], vec![(0, false)]);
        // the same bases on the plus strand read the node backwards
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let c = locus("c3", 1, '+', &[(1, 16)], &rc(S16), "C");
        let r = paf(&a, 0, 16, true, &c, 0, 16, "16M");
        let g = graph(vec![a, c], vec![r]);
        assert_eq!(g.paths[1], vec![(0, true)]);
    }

    #[test]
    fn an_internal_repeat_is_a_node_walked_twice() {
        let a = locus("c1", 1, '+', &[(1, 6), (11, 16)], "GATTACCCCCGATTAC", "A");
        let b = locus("c2", 1, '+', &[(1, 6)], "GATTAC", "B");
        let r1 = paf(&b, 0, 6, false, &a, 0, 6, "6M");
        let r2 = paf(&b, 0, 6, false, &a, 10, 16, "6M");
        let g = graph(vec![a, b], vec![r1, r2]);
        assert_eq!(g.nodes.len(), 1);
        assert_eq!(g.paths[0], vec![(0, false), (0, false)]);
        assert_eq!(g.stats.loops, 6);
    }

    #[test]
    fn same_dna_nodes_are_flagged() {
        let s = "ACGTACGTACGGATTACAGG";
        let a = locus("c1", 1, '+', &[(1, 20)], s, "A");
        let b = locus("c1", 11, '+', &[(1, 20)], &format!("{}TTTTTTTTTT", &s[10..]), "B");
        let r = paf(&a, 10, 20, false, &b, 0, 10, "10M");
        let g = graph(vec![a, b], vec![r]);
        let sd: Vec<&Node> = g.nodes.iter().filter(|n| n.same_dna).collect();
        assert_eq!(sd.len(), 1);
        assert_eq!((sd[0].seq.len(), sd[0].copies), (10, 2));
    }

    #[test]
    fn links_are_canonical_in_both_directions() {
        assert_eq!(canon(3, false, 1, false), canon(1, true, 3, true));
    }
```

- [ ] **Step 2: Run the unit tests and confirm they fail to compile**

Expected: `cannot find function build`, `cannot find type Graph`.

- [ ] **Step 3: Implement the compaction section**

Insert after the `// Joining` section:

```rust
// ------------------------------------------------------------------------------------------------------------------
// Compaction
// ------------------------------------------------------------------------------------------------------------------

/// A node: consecutive columns every path through them traverses together, cut at exon boundaries.
#[derive(Clone, Debug)]
struct Node {
    /// the representative copy's bases, in the node's orientation
    seq: Vec<u8>,
    /// distinct copies whose path traverses the node
    copies: usize,
    /// traversals, counting a copy that walks it twice (a loop) twice
    traversals: usize,
    /// mean and minimum identity against the node's sequence, over all traversals but the representative's own
    id_mean: f64,
    id_min: f64,
    /// every column is same-DNA (overlapping loci), not paralogy
    same_dna: bool,
}

/// A link between node ends, canonical so (a+, b+) and (b-, a-) are one link.
type Link = (usize, bool, usize, bool);

fn canon(a: usize, ar: bool, b: usize, br: bool) -> Link {
    let x = (a, ar, b, br);
    let y = (b, !br, a, !ar);
    if x <= y {
        x
    } else {
        y
    }
}

/// A graph view of one family.
#[derive(Clone, Debug)]
struct Graph {
    name: String,
    /// members sorted by key; `paths[i]` belongs to `members[i]`
    members: Vec<Locus>,
    /// node `n` is written as GFA id `n + 1`, numbered by first appearance along the paths
    nodes: Vec<Node>,
    /// per member: (node, reversed), in transcript order
    paths: Vec<Vec<(usize, bool)>>,
    links: BTreeSet<Link>,
    /// per member: bases differing from the node sequences it traverses
    mism: Vec<usize>,
    stats: Stats,
}

fn oriented(seq: &[u8], rev: bool) -> Vec<u8> {
    if rev {
        seq.iter().rev().map(|&b| comp(b)).collect()
    } else {
        seq.to_vec()
    }
}

fn build(name: &str, members: Vec<Locus>, pafs: &[&Paf]) -> Result<Graph> {
    let j = join(name, members, pafs)?;
    let ncol = j.rep.len();
    // column-level paths in transcript order: (column, reversed, exon index, vertex)
    let mut cpaths: Vec<Vec<(usize, u8, usize, u32)>> = Vec::new();
    for (i, m) in j.members.iter().enumerate() {
        let minus = (m.strand == '-') as u8;
        let mut steps: Vec<(usize, u8, usize, u32)> = j.exonic[i]
            .iter()
            .map(|&(o, e)| {
                let v = j.vid[i][o];
                (j.col[v as usize], j.orient[v as usize] ^ minus, e, v)
            })
            .collect();
        if minus == 1 {
            steps.reverse();
        }
        cpaths.push(steps);
    }
    // column ends: 2c = left end, 2c + 1 = right end, in the column's canonical orientation
    let exit = |c: usize, o: u8| 2 * c + if o == 0 { 1 } else { 0 };
    let entry = |c: usize, o: u8| 2 * c + if o == 0 { 0 } else { 1 };
    let mut nb: Vec<Vec<usize>> = vec![Vec::new(); 2 * ncol];
    let mut cut = vec![false; 2 * ncol];
    for p in &cpaths {
        if let (Some(f), Some(l)) = (p.first(), p.last()) {
            cut[entry(f.0, f.1)] = true;
            cut[exit(l.0, l.1)] = true;
        }
        for w in p.windows(2) {
            let (s, t) = (exit(w[0].0, w[0].1), entry(w[1].0, w[1].1));
            nb[s].push(t);
            nb[t].push(s);
            if w[0].2 != w[1].2 {
                cut[s] = true;
                cut[t] = true;
            }
        }
    }
    for v in nb.iter_mut() {
        v.sort_unstable();
        v.dedup();
    }
    // c's right end continues into d's left end when that is the only link on both ends and neither end is cut
    let next = |c: usize| -> Option<usize> {
        let s = 2 * c + 1;
        if cut[s] || nb[s].len() != 1 {
            return None;
        }
        let t = nb[s][0];
        if t % 2 != 0 || cut[t] || nb[t].len() != 1 || t / 2 == c {
            return None;
        }
        Some(t / 2)
    };
    let mut prev = vec![usize::MAX; ncol];
    let mut nxt = vec![usize::MAX; ncol];
    for c in 0..ncol {
        if let Some(d) = next(c) {
            nxt[c] = d;
            prev[d] = c;
        }
    }
    let mut node_of = vec![usize::MAX; ncol];
    let mut pos_in = vec![0usize; ncol];
    let mut chains: Vec<Vec<usize>> = Vec::new();
    for c in 0..ncol {
        if node_of[c] != usize::MAX {
            continue;
        }
        // walk back to the chain's head; a cycle of columns is broken where the walk returns to c
        let mut h = c;
        while prev[h] != usize::MAX && prev[h] != c && node_of[prev[h]] == usize::MAX {
            h = prev[h];
        }
        let mut chain = Vec::new();
        let mut x = h;
        loop {
            node_of[x] = chains.len();
            pos_in[x] = chain.len();
            chain.push(x);
            let y = nxt[x];
            if y == usize::MAX || node_of[y] != usize::MAX {
                break;
            }
            x = y;
        }
        chains.push(chain);
    }
    // a column's base in its canonical orientation: the representative's forward genomic base
    let base = |c: usize| -> u8 {
        let (m, o) = j.owner[j.rep[c] as usize];
        j.members[m].seq[o]
    };
    // the base a copy's transcript carries at vertex v
    let emitted = |i: usize, v: u32| -> u8 {
        let (_, o) = j.owner[v as usize];
        let b = j.members[i].seq[o];
        if j.members[i].strand == '-' {
            comp(b)
        } else {
            b
        }
    };
    // node-level paths; each traversal must enter a node at its first column (forward) or its last (reversed)
    let mut npaths: Vec<Vec<(usize, bool)>> = Vec::new();
    let mut mism = vec![0usize; j.members.len()];
    let mut trav_by_node: Vec<Vec<(usize, usize, usize, bool)>> = vec![Vec::new(); chains.len()]; // (member, len, mismatches, defining)
    for (i, p) in cpaths.iter().enumerate() {
        let who = j.members[i].path_name();
        let mut steps = Vec::new();
        let mut k = 0;
        while k < p.len() {
            let (c, o, _, _) = p[k];
            let nd = node_of[c];
            let chain = &chains[nd];
            let len = chain.len();
            let rev = o == 1;
            if pos_in[c] != if rev { len - 1 } else { 0 } {
                bail!("family {name}: path {who} enters node {} mid-way", nd + 1);
            }
            let (mut mm, mut defining) = (0usize, false);
            for x in 0..len {
                let &(c2, o2, _, v) =
                    p.get(k + x).with_context(|| format!("family {name}: path {who} ends inside node {}", nd + 1))?;
                let want = if rev { chain[len - 1 - x] } else { chain[x] };
                if c2 != want || (o2 == 1) != rev {
                    bail!("family {name}: path {who} leaves node {} mid-way", nd + 1);
                }
                let expect = if rev { comp(base(c2)) } else { base(c2) };
                if emitted(i, v) != expect {
                    mm += 1;
                }
                if v == j.rep[chain[0]] {
                    defining = true;
                }
            }
            if defining && mm > 0 {
                bail!("family {name}: the representative of node {} does not spell it (orientation bug)", nd + 1);
            }
            trav_by_node[nd].push((i, len, mm, defining));
            mism[i] += mm;
            steps.push((nd, rev));
            k += len;
        }
        npaths.push(steps);
    }
    // number nodes by first appearance along the paths (every column holds a vertex, so every node is on a path)
    let mut new_id = vec![usize::MAX; chains.len()];
    let mut order = Vec::new();
    for p in &npaths {
        for &(nd, _) in p {
            if new_id[nd] == usize::MAX {
                new_id[nd] = order.len();
                order.push(nd);
            }
        }
    }
    let paths: Vec<Vec<(usize, bool)>> =
        npaths.iter().map(|p| p.iter().map(|&(nd, r)| (new_id[nd], r)).collect()).collect();
    let mut nodes = Vec::with_capacity(order.len());
    for &nd in &order {
        let t = &trav_by_node[nd];
        let copies = t.iter().map(|x| x.0).collect::<BTreeSet<_>>().len();
        let ids: Vec<f64> = t.iter().filter(|x| !x.3).map(|x| 1.0 - x.2 as f64 / x.1 as f64).collect();
        let (id_mean, id_min) = if ids.is_empty() {
            (1.0, 1.0)
        } else {
            (ids.iter().sum::<f64>() / ids.len() as f64, ids.iter().cloned().fold(f64::INFINITY, f64::min))
        };
        nodes.push(Node {
            seq: chains[nd].iter().map(|&c| base(c)).collect(),
            copies,
            traversals: t.len(),
            id_mean,
            id_min,
            same_dna: chains[nd].iter().all(|&c| j.same_dna[c]),
        });
    }
    let mut links = BTreeSet::new();
    for p in &paths {
        for w in p.windows(2) {
            links.insert(canon(w[0].0, w[0].1, w[1].0, w[1].1));
        }
    }
    let stats = j.stats.clone();
    Ok(Graph { name: name.to_string(), members: j.members, nodes, paths, links, mism, stats })
}

/// The run-time invariants of the spec, checked on every graph before it is written.
fn check(g: &Graph) -> Result<()> {
    for (i, m) in g.members.iter().enumerate() {
        let who = m.path_name();
        let p = &g.paths[i];
        let len: usize = p.iter().map(|&(n, _)| g.nodes[n].seq.len()).sum();
        if len != m.exonic_bp() {
            bail!("family {}: path {who} spans {len} bp, its exonic length is {}", g.name, m.exonic_bp());
        }
        for w in p.windows(2) {
            if !g.links.contains(&canon(w[0].0, w[0].1, w[1].0, w[1].1)) {
                bail!("family {}: path {who} steps from node {} to {} without a link", g.name, w[0].0 + 1, w[1].0 + 1);
            }
        }
        let tx = m.transcript();
        let (mut diff, mut pos) = (0usize, 0usize);
        for &(n, r) in p {
            let s = oriented(&g.nodes[n].seq, r);
            let d = s.iter().zip(&tx[pos..pos + s.len()]).filter(|(a, b)| a != b).count();
            if d > 0 && g.nodes[n].traversals == 1 {
                bail!("family {}: path {who} differs from its transcript inside node {}, which only it walks", g.name, n + 1);
            }
            diff += d;
            pos += s.len();
        }
        if diff != g.mism[i] {
            bail!("family {}: path {who} spells {diff} differences, the identity tags count {}", g.name, g.mism[i]);
        }
    }
    let mut through: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); g.nodes.len()];
    for (i, p) in g.paths.iter().enumerate() {
        for &(n, _) in p {
            through[n].insert(i);
        }
    }
    for (n, node) in g.nodes.iter().enumerate() {
        if through[n].len() != node.copies {
            bail!("family {}: node {} has CN {} but {} paths walk it", g.name, n + 1, node.copies, through[n].len());
        }
    }
    Ok(())
}
```

- [ ] **Step 4: Run all tests**

Unit and integration commands. Expected: 8 new unit tests PASS, all earlier tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bin/family_gfa.rs
git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -F - <<'EOF'
family_gfa: compact columns into exon-aligned nodes, paths, tags, invariants (task 3/6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
EOF
```

---

### Task 4: Writers and the family run

**Files:**
- Modify: `src/bin/family_gfa.rs` (add a `// Writers` section; replace `main`; add unit tests)
- Modify: `tests/family_gfa.rs` (add end-to-end tests)

**Interfaces:**
- Consumes: `load`, `Inputs`, `build`, `check`, `Graph`, `Node`, `Link` from Tasks 1-3.
- Produces (used by Task 5): `const COMMENT: &str`; `fn node_name(prefix: &str, n: usize) -> String` (`n + 1`, prefixed
  `"{prefix}_"` when `prefix` is not empty); `fn sign(rev: bool) -> char`;
  `fn segments_and_links(out: &mut String, g: &Graph, prefix: &str)`;
  `fn path_line(name: &str, steps: &[(usize, bool)], prefix: &str) -> String`;
  `fn colour(copies: usize, n_copies: usize) -> String`;
  `fn colour_rows(g: &Graph, prefix: &str, red: &BTreeSet<String>) -> String` (rows only, no header);
  `fn family_graph(inp: &Inputs, family: &str) -> Result<Graph>` (built and checked);
  `fn write_file(path: &Path, text: &str) -> Result<()>`.

- [ ] **Step 1: Write the failing tests**

Add to `mod tests` in `src/bin/family_gfa.rs`:

```rust
    #[test]
    fn colours_follow_the_spec_ramp() {
        assert_eq!(colour(1, 5), "#eb6834");
        assert_eq!(colour(2, 5), "#86b6ef");
        assert_eq!(colour(5, 5), "#104281");
        assert_eq!(colour(2, 2), "#104281");
    }

    #[test]
    fn a_family_gfa_starts_with_the_header_and_the_comment() {
        let a = locus("c1", 1, '+', &[(1, 16)], S16, "A");
        let b = locus("c2", 1, '+', &[(1, 16)], S16, "B");
        let r = paf(&a, 0, 16, false, &b, 0, 16, "16M");
        let g = graph(vec![a, b], vec![r]);
        let text = family_gfa(&g);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "H\tVN:Z:1.0");
        assert_eq!(
            lines[1],
            "# graph view of O1 family F, built by family_gfa from O1's own alignments; not used to define families"
        );
        assert_eq!(lines[2], format!("S\t1\t{S16}\tLN:i:16\tCN:i:2\tID:f:1.0000\tIM:f:1.0000"));
        assert!(text.contains("P\tA@c1:1-16\t1+\t*\n"));
    }
```

Append to `tests/family_gfa.rs` (remove the `unused_helpers_kept_for_later_tasks` function from Task 1):

```rust
/// Every path step has an S line and every consecutive pair of steps has an L line (in either direction).
fn assert_gfa_consistent(text: &str) {
    let mut seg = HashSet::new();
    let mut links = HashSet::new();
    let mut paths = Vec::new();
    for l in text.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        match f[0] {
            "S" => {
                seg.insert(f[1].to_string());
            }
            "L" => {
                links.insert((f[1].to_string(), f[2].to_string(), f[3].to_string(), f[4].to_string()));
            }
            "P" => paths.push(f[2].to_string()),
            _ => {}
        }
    }
    let flip = |s: &str| if s == "+" { "-".to_string() } else { "+".to_string() };
    for p in paths {
        let steps: Vec<(String, String)> =
            p.split(',').map(|s| (s[..s.len() - 1].to_string(), s[s.len() - 1..].to_string())).collect();
        for (n, _) in &steps {
            assert!(seg.contains(n), "path step {n} has no S line");
        }
        for w in steps.windows(2) {
            let fwd = (w[0].0.clone(), w[0].1.clone(), w[1].0.clone(), w[1].1.clone());
            let bwd = (w[1].0.clone(), flip(&w[1].1), w[0].0.clone(), flip(&w[0].1));
            assert!(links.contains(&fwd) || links.contains(&bwd), "no L line for {fwd:?}");
        }
    }
}

/// `vg` parses the GFA; None when vg is not installed.
fn vg_parses(path: &Path) -> Option<bool> {
    let o = Command::new("vg").args(["convert", "-g", path.to_str().unwrap(), "-p"]).output().ok()?;
    Some(o.status.success())
}

fn read(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

#[test]
fn writes_one_graph_per_family_with_the_expected_structure() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let out = d.path().join("out");
    let o = run(&["--prefix", &p, "--outdir", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", stderr(&o));
    let f1 = read(&out.join("FAM1.gfa"));
    assert!(f1.starts_with("H\tVN:Z:1.0\n# graph view of O1 family FAM1, "));
    // two exon nodes shared by all three copies; B's SNP sits in exon 1
    assert!(f1.contains("\tLN:i:20\tCN:i:3\tID:f:0.9750\tIM:f:0.9500\n"), "{f1}");
    assert!(f1.contains("\tLN:i:20\tCN:i:3\tID:f:1.0000\tIM:f:1.0000\n"), "{f1}");
    // C is the minus-strand reverse complement of A: its transcript walks the same nodes forward
    assert!(f1.contains("P\tGA@chr1:1001-1060\t1+,2+\t*\n"), "{f1}");
    assert!(f1.contains("P\tGC@chr3:3001-3060\t1+,2+\t*\n"), "{f1}");
    assert_gfa_consistent(&f1);
    assert_gfa_consistent(&read(&out.join("FAM2.gfa")));
    let summary = read(&out.join("summary.tsv"));
    assert!(summary.starts_with(
        "family\tn_copies\tn_nodes\tn_shared_nodes\tshared_exonic_frac\tloops\tinversion_conflicts\tsame_dna_bp\n"
    ));
    assert!(summary.contains("FAM1\t3\t2\t2\t1.0000\t0\t0\t0\n"), "{summary}");
    assert!(summary.contains("FAM2\t2\t1\t1\t1.0000\t0\t0\t0\n"), "{summary}");
    let paths = read(&out.join("FAM1.paths.tsv"));
    assert!(paths.contains("GC@chr3:3001-3060\tGC\tchr3:3001-3060\t-\t2\t40\n"), "{paths}");
    let colours = read(&out.join("FAM1.colours.csv"));
    assert!(colours.starts_with("Node,Colour\n1,#104281\n"), "{colours}");
    match vg_parses(&out.join("FAM1.gfa")) {
        Some(ok) => assert!(ok, "vg could not parse FAM1.gfa"),
        None => eprintln!("vg not installed: GFA parse check skipped"),
    }
}

#[test]
fn reruns_are_byte_identical() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let (o1, o2) = (d.path().join("o1"), d.path().join("o2"));
    assert!(run(&["--prefix", &p, "--outdir", o1.to_str().unwrap()]).status.success());
    assert!(run(&["--prefix", &p, "--outdir", o2.to_str().unwrap()]).status.success());
    for f in ["FAM1.gfa", "FAM1.colours.csv", "FAM1.paths.tsv", "FAM2.gfa", "summary.tsv"] {
        assert_eq!(read(&o1.join(f)), read(&o2.join(f)), "{f} differs between reruns");
    }
}

#[test]
fn families_limits_the_output() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let out = d.path().join("out");
    assert!(run(&["--prefix", &p, "--outdir", out.to_str().unwrap(), "--families", "FAM2"]).status.success());
    assert!(out.join("FAM2.gfa").exists());
    assert!(!out.join("FAM1.gfa").exists());
}
```

Also change the Task 1 test `loads_a_valid_run` assertion to
`assert!(stderr(&o).contains("2 families written"), "{}", stderr(&o));`.

- [ ] **Step 2: Run the tests and confirm they fail**

Expected: compile error `cannot find function colour` / `family_gfa`; integration tests fail (no output files).

- [ ] **Step 3: Implement the writers and the new `main`**

Insert after the `// Compaction` section:

```rust
// ------------------------------------------------------------------------------------------------------------------
// Writers
// ------------------------------------------------------------------------------------------------------------------

const COMMENT: &str = "built by family_gfa from O1's own alignments; not used to define families";

fn node_name(prefix: &str, n: usize) -> String {
    if prefix.is_empty() {
        (n + 1).to_string()
    } else {
        format!("{prefix}_{}", n + 1)
    }
}

fn sign(rev: bool) -> char {
    if rev {
        '-'
    } else {
        '+'
    }
}

/// The S and L lines of one graph; node names bare (family GFA) or prefixed (joint GFA).
fn segments_and_links(out: &mut String, g: &Graph, prefix: &str) {
    for (n, node) in g.nodes.iter().enumerate() {
        out.push_str(&format!(
            "S\t{}\t{}\tLN:i:{}\tCN:i:{}\tID:f:{:.4}\tIM:f:{:.4}",
            node_name(prefix, n),
            String::from_utf8_lossy(&node.seq),
            node.seq.len(),
            node.copies,
            node.id_mean,
            node.id_min
        ));
        if node.same_dna {
            out.push_str("\tSD:i:1");
        }
        out.push('\n');
    }
    for &(a, ar, b, br) in &g.links {
        out.push_str(&format!("L\t{}\t{}\t{}\t{}\t0M\n", node_name(prefix, a), sign(ar), node_name(prefix, b), sign(br)));
    }
}

fn path_line(name: &str, steps: &[(usize, bool)], prefix: &str) -> String {
    let s: Vec<String> = steps.iter().map(|&(n, r)| format!("{}{}", node_name(prefix, n), sign(r))).collect();
    format!("P\t{name}\t{}\t*\n", s.join(","))
}

fn family_gfa(g: &Graph) -> String {
    let mut out = format!("H\tVN:Z:1.0\n# graph view of O1 family {}, {COMMENT}\n", g.name);
    segments_and_links(&mut out, g, "");
    for (i, m) in g.members.iter().enumerate() {
        out.push_str(&path_line(&m.path_name(), &g.paths[i], ""));
    }
    out
}

/// Bandage colour: private nodes orange; shared nodes light blue (2 copies) to dark blue (all copies).
fn colour(copies: usize, n_copies: usize) -> String {
    if copies <= 1 {
        return "#eb6834".to_string();
    }
    let t = if n_copies <= 2 { 1.0 } else { (copies - 2) as f64 / (n_copies - 2) as f64 };
    let (lo, hi) = ([0x86u8, 0xb6, 0xef], [0x10u8, 0x42, 0x81]);
    let mix = |i: usize| (lo[i] as f64 + (hi[i] as f64 - lo[i] as f64) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", mix(0), mix(1), mix(2))
}

/// Colour rows (no header) for one graph; nodes named in `red` are fusion-link ends.
fn colour_rows(g: &Graph, prefix: &str, red: &BTreeSet<String>) -> String {
    let mut out = String::new();
    for (n, node) in g.nodes.iter().enumerate() {
        let name = node_name(prefix, n);
        let c = if red.contains(&name) { "#e34948".to_string() } else { colour(node.copies, g.members.len()) };
        out.push_str(&format!("{name},{c}\n"));
    }
    out
}

fn paths_tsv(g: &Graph) -> String {
    let mut out = String::from("path_name\tgene_id\tlocus\tstrand\tn_exon\texonic_bp\n");
    for m in &g.members {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            m.path_name(),
            m.gene_id,
            m.key.name(),
            m.strand,
            m.exons.len(),
            m.exonic_bp()
        ));
    }
    out
}

const SUMMARY_HEADER: &str =
    "family\tn_copies\tn_nodes\tn_shared_nodes\tshared_exonic_frac\tloops\tinversion_conflicts\tsame_dna_bp\n";

fn summary_row(g: &Graph) -> String {
    let total: usize = g.members.iter().map(|m| m.exonic_bp()).sum();
    let shared: usize = g
        .paths
        .iter()
        .flat_map(|p| p.iter())
        .filter(|&&(n, _)| g.nodes[n].copies >= 2)
        .map(|&(n, _)| g.nodes[n].seq.len())
        .sum();
    format!(
        "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\n",
        g.name,
        g.members.len(),
        g.nodes.len(),
        g.nodes.iter().filter(|n| n.copies >= 2).count(),
        if total == 0 { 0.0 } else { shared as f64 / total as f64 },
        g.stats.loops,
        g.stats.inversion_conflicts,
        g.stats.same_dna_bp
    )
}

fn write_file(path: &Path, text: &str) -> Result<()> {
    std::fs::write(path, text).with_context(|| format!("writing {}", path.display()))
}

/// The family's graph, built from its members and the records between them, and checked.
fn family_graph(inp: &Inputs, family: &str) -> Result<Graph> {
    let keys = &inp.all[family];
    let members: Vec<Locus> = keys.iter().map(|k| inp.loci[k].clone()).collect();
    let set: HashSet<&Key> = keys.iter().collect();
    let pafs: Vec<&Paf> = inp.pafs.iter().filter(|r| set.contains(&r.q) && set.contains(&r.t)).collect();
    let g = build(family, members, &pafs)?;
    check(&g)?;
    Ok(g)
}
```

Replace `main` with:

```rust
fn main() -> Result<()> {
    let a = Args::parse();
    let inp = load(&a.prefix, &a.families, &[])?;
    std::fs::create_dir_all(&a.outdir).with_context(|| format!("creating {}", a.outdir.display()))?;
    let mut summary = String::from(SUMMARY_HEADER);
    for f in &inp.selected {
        let g = family_graph(&inp, f)?;
        write_file(&a.outdir.join(format!("{f}.gfa")), &family_gfa(&g))?;
        write_file(
            &a.outdir.join(format!("{f}.colours.csv")),
            &format!("Node,Colour\n{}", colour_rows(&g, "", &BTreeSet::new())),
        )?;
        write_file(&a.outdir.join(format!("{f}.paths.tsv")), &paths_tsv(&g))?;
        summary.push_str(&summary_row(&g));
    }
    write_file(&a.outdir.join("summary.tsv"), &summary)?;
    eprintln!(
        "family_gfa: {} families written to {} ({} alignment records, {} without cg:Z skipped)",
        inp.selected.len(),
        a.outdir.display(),
        inp.pafs.len(),
        inp.no_cigar
    );
    Ok(())
}
```

- [ ] **Step 4: Run all tests**

Expected: all unit and integration tests PASS. If `vg` is installed (`/home/juanfra/miniforge3/bin/vg`), the vg
check runs; it prints a warning about the `#` comment line and must still exit 0.

- [ ] **Step 5: Commit**

```bash
git add src/bin/family_gfa.rs tests/family_gfa.rs
git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -F - <<'EOF'
family_gfa: GFA, colours, paths and summary writers; family run end to end (task 4/6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
EOF
```

---

### Task 5: Fusions across families

**Files:**
- Modify: `src/bin/family_gfa.rs` (add a `// Fusions` section, a `--relations` argument, fusion output in `main`)
- Modify: `tests/family_gfa.rs` (fusion tests)

**Interfaces:**
- Consumes: `load` (with `extra` keys), `family_graph`, `build`, `check`, `segments_and_links`, `path_line`,
  `node_name`, `sign`, `colour_rows`, `write_file`, `COMMENT` from Tasks 1-4.
- Produces: `struct Relation { id: String, transcript: String, units: Vec<Key>, listed: Vec<String> }`;
  `fn read_relations(path: &Path) -> Result<Vec<Relation>>`;
  `struct FusionOut { families: Vec<String>, note: String, reason: String, gfa: Option<(String, String)> }`;
  `fn fusion(inp: &Inputs, r: &Relation, cache: &mut BTreeMap<String, Graph>) -> Result<FusionOut>`.

- [ ] **Step 1: Write the failing integration tests**

Append to `tests/family_gfa.rs`:

```rust
const REL_HEADER: &str = "relation_id\ttranscript\tfused_gene_id\tfused_locus\tstrand\tn_units\tcut_introns\tdetector\tunit_ids\tunit_loci\tunit_families\tunit_family_sizes\toutcome\trelation\tseparated\tref_lenient\tref_strict\n";

fn rel_row(id: &str, tx: &str, units: &[&str], fams: &[&str]) -> String {
    format!(
        "{id}\t{tx}\t{tx}\t.\t+\t{}\t.\tf1\t.\t{}\t{}\t.\t.\t.\ttrue\t.\t.\n",
        units.len(),
        units.join(","),
        fams.join(",")
    )
}

#[test]
fn fusions_get_joint_graphs_and_reasons() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let (a, b, dd, f) = (l[0].key(), l[1].key(), l[3].key(), l[5].key());
    let rels = format!(
        "{REL_HEADER}{}{}{}{}{}",
        rel_row("REL1", "TX1", &[&a, &dd], &["FAM1", "FAM2"]),
        rel_row("REL2", "TX2", &[&a, &b], &["FAM1", "FAM1"]),
        rel_row("REL3", "TX3", &[&a, &f, &dd], &["FAM1", "-", "FAM2"]),
        rel_row("REL4", "TX4", &["chr9:1-10", &dd], &["-", "FAM2"]),
        rel_row("REL5", "TX5", &[&a, &dd], &["FAM2", "FAM2"]),
    );
    let rp = d.path().join("t.relations.tsv");
    std::fs::write(&rp, rels).unwrap();
    let out = d.path().join("out");
    let o = run(&["--prefix", &p, "--outdir", out.to_str().unwrap(), "--relations", rp.to_str().unwrap()]);
    assert!(o.status.success(), "{}", stderr(&o));

    let j1 = read(&out.join("fusion_REL1.gfa"));
    assert!(j1.starts_with("H\tVN:Z:1.0\n# graph view of fusion REL1 (TX1) across O1 families FAM1, FAM2, "), "{j1}");
    assert!(j1.contains("P\tTX1|fused\tFAM1_1+,FAM1_2+,FAM2_1+\t*\n"), "{j1}");
    assert!(j1.contains("L\tFAM1_2\t+\tFAM2_1\t+\t0M\tFU:Z:REL1\n"), "{j1}");
    assert!(j1.contains("P\tFAM1|GA@chr1:1001-1060\tFAM1_1+,FAM1_2+\t*\n"), "{j1}");
    assert_gfa_consistent(&j1);
    let c1 = read(&out.join("fusion_REL1.colours.csv"));
    assert!(c1.contains("FAM1_2,#e34948\n") && c1.contains("FAM2_1,#e34948\n"), "{c1}");

    let j3 = read(&out.join("fusion_REL3.gfa"));
    assert!(j3.contains("P\tTX3|fused\tFAM1_1+,FAM1_2+,U2_1+,FAM2_1+\t*\n"), "{j3}");
    assert_gfa_consistent(&j3);

    let fus = read(&out.join("fusions.tsv"));
    assert!(fus.starts_with("relation_id\ttranscript\tfamilies\tunits\tjoint_gfa\treason\tnote\n"), "{fus}");
    assert!(fus.contains("REL1\tTX1\tFAM1,FAM2\t"), "{fus}");
    assert!(fus.lines().any(|x| x.starts_with("REL2\t") && x.contains("\t-\tall units are in one family\t")), "{fus}");
    assert!(fus.lines().any(|x| x.starts_with("REL4\t") && x.contains("unit chr9:1-10 is not in the loci files")), "{fus}");
    let rel5 = fus.lines().find(|x| x.starts_with("REL5\t")).unwrap();
    assert!(rel5.contains("fusion_REL5.gfa") && rel5.contains("chr1:1001-1060 listed FAM2, clusters.tsv FAM1"), "{rel5}");
    assert!(!out.join("fusion_REL2.gfa").exists());
}
```

- [ ] **Step 2: Run the integration tests and confirm the new one fails**

Expected: FAIL: clap rejects `--relations` (unexpected argument).

- [ ] **Step 3: Implement the fusions section**

Add the argument to `Args` (after `families`):

```rust
    /// `mcl_families --emit-relations` output: writes fusions.tsv, and a joint GFA per fusion across families.
    #[arg(long)]
    relations: Option<PathBuf>,
```

Insert after the `// Writers` section:

```rust
// ------------------------------------------------------------------------------------------------------------------
// Fusions
// ------------------------------------------------------------------------------------------------------------------

/// A fused transcript from `mcl_families --emit-relations`, its units in transcript order, and the family each unit
/// was listed with ("-" = none).
struct Relation {
    id: String,
    transcript: String,
    units: Vec<Key>,
    listed: Vec<String>,
}

fn read_relations(path: &Path) -> Result<Vec<Relation>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().with_context(|| format!("{} is empty", path.display()))??;
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| {
        cols.iter().position(|c| *c == name).with_context(|| format!("{} has no {name} column", path.display()))
    };
    let (ii, ti, ui, fi) = (idx("relation_id")?, idx("transcript")?, idx("unit_loci")?, idx("unit_families")?);
    let mut out = Vec::new();
    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let units: Vec<Key> = field(&f, ui, &line)?.split(',').map(Key::parse).collect::<Result<_>>()?;
        let listed: Vec<String> = field(&f, fi, &line)?.split(',').map(|s| s.to_string()).collect();
        if units.len() != listed.len() {
            bail!("{}: unit_loci and unit_families differ in length: {line}", path.display());
        }
        out.push(Relation {
            id: field(&f, ii, &line)?.to_string(),
            transcript: field(&f, ti, &line)?.to_string(),
            units,
            listed,
        });
    }
    Ok(out)
}

/// One fusions.tsv row, plus the joint GFA and its colours when the relation spans two or more families.
struct FusionOut {
    families: Vec<String>,
    note: String,
    reason: String,
    gfa: Option<(String, String)>,
}

/// Units are assigned to families by `.clusters.tsv` (a disagreeing `unit_families` entry goes to `note`). A relation
/// whose units span two or more families gets the joint graph: each family's graph (node names prefixed by the family),
/// a unit in no family as its own graph `U<n>`, and the fused transcript as a path joining its units' paths by
/// `FU:Z:`-tagged links.
fn fusion(inp: &Inputs, r: &Relation, cache: &mut BTreeMap<String, Graph>) -> Result<FusionOut> {
    let fams: Vec<Option<String>> = r.units.iter().map(|k| inp.member_of.get(k).cloned()).collect();
    let names: Vec<String> = fams.iter().map(|f| f.clone().unwrap_or_else(|| "-".into())).collect();
    let note: Vec<String> = r
        .units
        .iter()
        .zip(names.iter().zip(&r.listed))
        .filter(|(_, (got, listed))| got != listed)
        .map(|(k, (got, listed))| format!("{} listed {listed}, clusters.tsv {got}", k.name()))
        .collect();
    let note = note.join("; ");
    let none = |reason: String| FusionOut { families: names.clone(), note: note.clone(), reason, gfa: None };
    if let Some(k) = r.units.iter().find(|k| !inp.loci.contains_key(*k)) {
        return Ok(none(format!("unit {} is not in the loci files", k.name())));
    }
    let distinct: BTreeSet<String> = fams.iter().flatten().cloned().collect();
    if distinct.is_empty() {
        return Ok(none("no unit is in a family".into()));
    }
    if distinct.len() == 1 {
        return Ok(none("all units are in one family".into()));
    }
    for f in &distinct {
        if !cache.contains_key(f) {
            let g = family_graph(inp, f)?;
            cache.insert(f.clone(), g);
        }
    }
    let mut unit_graphs: BTreeMap<String, Graph> = BTreeMap::new();
    for (u, (k, f)) in r.units.iter().zip(&fams).enumerate() {
        if f.is_none() {
            let name = format!("U{}", u + 1);
            let g = build(&name, vec![inp.loci[k].clone()], &[])?;
            check(&g)?;
            unit_graphs.insert(name, g);
        }
    }
    let mut graphs: Vec<(String, &Graph)> = distinct.iter().map(|f| (f.clone(), &cache[f])).collect();
    graphs.extend(unit_graphs.iter().map(|(n, g)| (n.clone(), g)));
    // the fused path: each unit's own path, in transcript order, joined by fusion links
    let mut fused: Vec<String> = Vec::new();
    let mut fusion_links = String::new();
    let mut red: BTreeSet<String> = BTreeSet::new();
    let mut prev_end: Option<(String, bool)> = None;
    for (u, (k, f)) in r.units.iter().zip(&fams).enumerate() {
        let prefix = f.clone().unwrap_or_else(|| format!("U{}", u + 1));
        let g = match f {
            Some(f) => &cache[f],
            None => &unit_graphs[&prefix],
        };
        let mi = g
            .members
            .iter()
            .position(|m| &m.key == k)
            .with_context(|| format!("relation {}: unit {} is not in its graph", r.id, k.name()))?;
        let steps = &g.paths[mi];
        let (first, last) = (steps[0], *steps.last().expect("a path has at least one exonic base"));
        if let Some((pn, pr)) = prev_end.take() {
            let next = node_name(&prefix, first.0);
            fusion_links.push_str(&format!(
                "L\t{pn}\t{}\t{next}\t{}\t0M\tFU:Z:{}\n",
                sign(pr),
                sign(first.1),
                r.id
            ));
            red.insert(pn);
            red.insert(next);
        }
        prev_end = Some((node_name(&prefix, last.0), last.1));
        fused.extend(steps.iter().map(|&(n, rv)| format!("{}{}", node_name(&prefix, n), sign(rv))));
    }
    let fam_list: Vec<String> = distinct.iter().cloned().collect();
    let mut text = format!(
        "H\tVN:Z:1.0\n# graph view of fusion {} ({}) across O1 families {}, {COMMENT}\n",
        r.id,
        r.transcript,
        fam_list.join(", ")
    );
    let mut colours = String::from("Node,Colour\n");
    let mut paths = String::new();
    for (prefix, g) in &graphs {
        segments_and_links(&mut text, g, prefix);
        colours.push_str(&colour_rows(g, prefix, &red));
        for (i, m) in g.members.iter().enumerate() {
            paths.push_str(&path_line(&format!("{prefix}|{}", m.path_name()), &g.paths[i], prefix));
        }
    }
    text.push_str(&fusion_links);
    text.push_str(&paths);
    text.push_str(&format!("P\t{}|fused\t{}\t*\n", r.transcript, fused.join(",")));
    Ok(FusionOut { families: names, note, reason: String::new(), gfa: Some((text, colours)) })
}
```

In `main`, replace `let inp = load(&a.prefix, &a.families, &[])?;` with:

```rust
    let rels = match &a.relations {
        Some(p) => read_relations(p)?,
        None => Vec::new(),
    };
    let extra: Vec<Key> = rels.iter().flat_map(|r| r.units.iter().cloned()).collect();
    let inp = load(&a.prefix, &a.families, &extra)?;
```

and, just before `eprintln!` in `main`, add:

```rust
    if a.relations.is_some() {
        let mut cache: BTreeMap<String, Graph> = BTreeMap::new();
        let mut tsv = String::from("relation_id\ttranscript\tfamilies\tunits\tjoint_gfa\treason\tnote\n");
        for r in &rels {
            let fo = fusion(&inp, r, &mut cache)?;
            let file = match &fo.gfa {
                Some((text, colours)) => {
                    let name = format!("fusion_{}.gfa", r.id);
                    write_file(&a.outdir.join(&name), text)?;
                    write_file(&a.outdir.join(format!("fusion_{}.colours.csv", r.id)), colours)?;
                    name
                }
                None => "-".to_string(),
            };
            let units: Vec<String> = r.units.iter().map(|k| k.name()).collect();
            tsv.push_str(&format!(
                "{}\t{}\t{}\t{}\t{file}\t{}\t{}\n",
                r.id,
                r.transcript,
                fo.families.join(","),
                units.join(","),
                fo.reason,
                fo.note
            ));
        }
        write_file(&a.outdir.join("fusions.tsv"), &tsv)?;
    }
```

- [ ] **Step 4: Run all tests**

Expected: all unit and integration tests PASS (including `fusions_get_joint_graphs_and_reasons`).

- [ ] **Step 5: Commit**

```bash
git add src/bin/family_gfa.rs tests/family_gfa.rs
git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -F - <<'EOF'
family_gfa: joint graphs for fusions across families, fusions.tsv (task 5/6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
EOF
```

---

### Task 6: Real-data smoke test and report

**Files:**
- Create: `docs/FAMILY_GFA_SMOKE_2026-09-30.md`

**Interfaces:**
- Consumes: the built `family_gfa` binary.
- Produces: the smoke report (numbers per flagship family and fusion, invariants passed, vg parse result).

Input run (read-only, written by another session at 2026-09-30 15:50 with current defaults, `--emit-units` and
`--emit-relations`): `P=/mnt/linuxdisk/tmp/rustle_figures_dev/units_port/H/f1units/hsa16.fam` (170 relations, 30
spanning two or more families; all 349 unit keys are loci in `$P.loci.gff3`). The spec's step "rerun the chr16
families stage" is met by this run, because `--emit-relations` exists only in that session's uncommitted code.

- [ ] **Step 1: Build the release binary**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2
export CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2
tools/rlock.sh heavy cargo build --release --bin family_gfa > /mnt/linuxdisk/tmp/fg_build.log 2>&1; echo exit=$?
```
Expected: `exit=0`; binary at `$CARGO_TARGET_DIR/release/family_gfa`.

- [ ] **Step 2: Name the families by RefSeq genes and pick the flagships**

```bash
P=/mnt/linuxdisk/tmp/rustle_figures_dev/units_port/H/f1units/hsa16.fam
tools/rlock.sh light python3 - "$P" > /mnt/linuxdisk/tmp/fg_names.tsv <<'PY'
import csv, gzip, sys
from collections import defaultdict
P = sys.argv[1]
genes = []
with gzip.open('/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz', 'rt') as fh:
    for l in fh:
        f = l.split('\t')
        if len(f) > 8 and f[0] == 'chr16' and f[2] == 'gene':
            name = [x[5:] for x in f[8].split(';') if x.startswith('Name=')]
            if name:
                genes.append((int(f[3]), int(f[4]), name[0]))
fams = defaultdict(list)
for r in csv.DictReader(open(P + '.clusters.tsv'), delimiter='\t'):
    s, e = int(r['start']), int(r['end'])
    hit = sorted({n for a, b, n in genes if a <= e and b >= s})
    fams[r['cluster_id']].append(','.join(hit) or '-')
print('family\tn\tgenes')
for f, ms in sorted(fams.items(), key=lambda x: (-len(x[1]), x[0])):
    print(f"{f}\t{len(ms)}\t{';'.join(ms)}")
npip = max(fams, key=lambda f: (sum('NPIP' in g for g in fams[f]), -int(f[3:]) if f[3:].isdigit() else 0))
pkd = sorted(f for f in fams if any('PKD1P' in g for g in fams[f]))
small = next((f for f, ms in sorted(fams.items()) if len(ms) == 2 and '-' not in ms and ms[0] != ms[1]), None)
chosen = []
for f in [npip] + pkd + ([small] if small else []):
    if f not in chosen:
        chosen.append(f)
rels = [r for r in csv.DictReader(open(P + '.relations.tsv'), delimiter='\t')
        if npip in r['unit_families'].split(',') and len(set(r['unit_families'].split(',')) - {'-'}) >= 2]
print('#FLAGSHIPS\t' + ','.join(chosen))
print(f"#NPIP\t{npip}\t#PKD1P\t{','.join(pkd) or '-'}\t#SMALL\t{small}")
print('#NPIP_CROSS_FAMILY_RELATIONS\t' + (','.join(r['relation_id'] for r in rels) or '-'))
PY
grep '^#' /mnt/linuxdisk/tmp/fg_names.tsv
FLAG=$(grep '^#FLAGSHIPS' /mnt/linuxdisk/tmp/fg_names.tsv | cut -f2); echo "flagships: $FLAG"
```
Expected: `#NPIP` names a family holding most NPIP* genes; `#PKD1P` lists the family or families holding PKD1P*
genes; `#NPIP_CROSS_FAMILY_RELATIONS` lists the relations joining NPIP to another family (the PKD1P6-NPIPP1 fusion is
among them if the fused transcript was cut into units). Copy all three lines into the report.

- [ ] **Step 3: Run the exporter on the flagships and all relations**

```bash
OUT=/mnt/linuxdisk/tmp/family_gfa_smoke
BIN=$CARGO_TARGET_DIR/release/family_gfa
/usr/bin/time -v $BIN --prefix $P --outdir $OUT --families "$FLAG" --relations $P.relations.tsv \
  > $OUT.log 2>&1; echo exit=$?
tail -3 $OUT.log; grep -E "Maximum resident|Elapsed" $OUT.log
cat $OUT/summary.tsv
awk -F'\t' 'NR>1 && $5!="-"' $OUT/fusions.tsv | wc -l
awk -F'\t' 'NR>1 {print $6}' $OUT/fusions.tsv | sort | uniq -c
```
Expected: `exit=0` (every family and joint graph passed the run-time invariants; otherwise the run stops with the
invariant's message); peak memory and time recorded; 30 joint graphs expected (the relations spanning two or more
families), give or take the 6 units whose listed family disagrees with `.clusters.tsv`.

- [ ] **Step 3b: TBC1D3, if a chr17 units run exists**

```bash
ls /mnt/linuxdisk/tmp/rustle_figures_dev/units_port/*/*/hsa17.fam.clusters.tsv 2>/dev/null || echo "no chr17 units run"
```
If a file is listed, repeat Steps 2-3 with `P` set to that prefix, `chr16` replaced by `chr17` in the script, and
`TBC1D3` in place of `NPIP` for the flagship pick; otherwise write "no chr17 units run on disk; TBC1D3 deferred" in the
report.

- [ ] **Step 4: Validate every GFA with vg and run the whole catalog once**

```bash
for g in $OUT/*.gfa; do vg convert -g "$g" -p > /dev/null 2>/mnt/linuxdisk/tmp/fg_vg.err || echo "vg FAILED $g"; done; echo vg-done
/usr/bin/time -v $BIN --prefix $P --outdir ${OUT}_all > ${OUT}_all.log 2>&1; echo exit=$?
grep -E "Maximum resident|Elapsed" ${OUT}_all.log; awk 'NR>1' ${OUT}_all/summary.tsv | wc -l
```
Expected: no `vg FAILED` lines; the all-families run exits 0 and writes one summary row per family.

- [ ] **Step 5: Write the report**

Create `docs/FAMILY_GFA_SMOKE_2026-09-30.md` with: the input run and why it was used (Task 6 header), the build
command, a table per flagship (`family, genes, n_copies, n_nodes, n_shared_nodes, shared_exonic_frac, loops,
inversion_conflicts, same_dna_bp` copied from `summary.tsv`), the PKD1P6-NPIPP1 relation id(s) and their joint GFA
file names, the count of joint graphs written and the `reason` counts from `fusions.tsv`, vg results, whole-catalog
time and peak memory, and how to view: "load `<family>.gfa` in Bandage, then File > Load CSV data >
`<family>.colours.csv`". State that the graphs are a graph view of O1 families, built from O1's own alignments, and not
used to define families.

- [ ] **Step 6: Commit**

```bash
git add docs/FAMILY_GFA_SMOKE_2026-09-30.md
git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -F - <<'EOF'
family_gfa: chr16 smoke test (NPIP, PKD1P6-NPIPP1 fusion, whole catalog) (task 6/6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb
EOF
```
