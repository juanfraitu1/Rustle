//! On-disk cache of the catalog builder's expensive intermediates, for fast re-runs and for analysis.
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `RUSTLE_CACHE_DIR`, set by default by `tools/rustle_pipeline.sh`)
//!
//! Enabled by `RUSTLE_CACHE_DIR=<dir>` (the pipeline driver sets `PREFIX.cache`); unset = nothing is read or
//! written and every output is byte-identical to a build without this module.
//!
//! Two objects are cached, both at boundaries where everything downstream reads only the cached object:
//!
//! * **`reps/<key>/`** — the collapsed locus representatives with their read statistics, i.e. the state of
//!   `detect_homology_catalog_genome_wide` after both BAM passes and the locus collapse (human chr16: ~260 of
//!   ~370 s). Files: `reps.tsv` (one row per representative, in representative-index order), `reps.fa` (the
//!   exact sequence bytes, one line each), `key.tsv` (the full key material), `DONE`.
//! * **`paf/<key>/`** — one all-vs-all minimap2 PAF per (query bytes, command line, minimap2 version):
//!   `out.paf`, `key.tsv`, `DONE`.
//!
//! A hit requires `DONE` and a `key.tsv` byte-identical to the key the current run computes, so a hash
//! collision in the directory name can only cause a miss, never a wrong hit. Writes go to a temporary
//! directory renamed into place, so an interrupted run leaves no partial entry.
//!
//! The representatives key covers: the executable itself (path, size, mtime — any rebuild invalidates), the
//! BAM and its index, the FASTA and its index (path, size, mtime), the `DenovoConfig`, and every `RUSTLE_*`
//! variable EXCEPT [`DOWNSTREAM_ONLY_ENV`], the settings read only after the boundary (each verified by grep
//! on 2026-09-24: the E_r edge rule, gamma, the coverage split, the edge dump and logging switches). An
//! unknown variable therefore over-invalidates; it can never produce a stale hit.
use crate::vg_family::family_detect::DenovoTranscript;
use anyhow::{Context, Result};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// `RUSTLE_*` settings read only downstream of the representatives boundary (or output-neutral there). Diagnostics
/// that print UPSTREAM (`RUSTLE_COLLAPSE_STATS`, `RUSTLE_LOCUS_AUDIT`, `RUSTLE_DEBUG_LOCUS`) are deliberately NOT
/// listed: setting one changes the key, so the representatives are recomputed and the diagnostics print. Varying them
/// re-uses the cached representatives — which is the point: edge-rule sweeps then skip both BAM passes and
/// the collapse. Anything else in the environment is part of the key.
pub const DOWNSTREAM_ONLY_ENV: &[&str] = &[
    "RUSTLE_CACHE_DIR",
    "RUSTLE_ER_MIN_COVERAGE",
    "RUSTLE_ER_COVERAGE_LONGER_FLOOR",
    "RUSTLE_ER_SENSITIVE_ONLY",
    "RUSTLE_ER_EDGE_DUMP",
    "RUSTLE_GENOME_GAMMA",
    "RUSTLE_COVERAGE_SPLIT",
    "RUSTLE_POA_MEMO",
];

/// The cache root, or `None` when caching is off.
pub fn cache_root() -> Option<PathBuf> {
    std::env::var("RUSTLE_CACHE_DIR").ok().filter(|v| !v.is_empty()).map(PathBuf::from)
}

/// FNV-1a 64 — stable across Rust releases (unlike `DefaultHasher`), used only to NAME entries; the key text
/// itself is compared in full on every hit.
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Incremental FNV-1a 64 (same values as [`fnv1a64`] over the concatenated input).
#[derive(Clone, Copy)]
pub struct Fnv(u64);
impl Default for Fnv {
    fn default() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
}
impl Fnv {
    pub fn update(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 ^= b as u64;
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    pub fn finish(self) -> u64 {
        self.0
    }
}

/// `path<TAB>size<TAB>mtime_ns` of a file (canonical path), or `path<TAB>absent`.
pub fn file_fingerprint(path: &str) -> String {
    let canon = std::fs::canonicalize(path).map(|p| p.display().to_string()).unwrap_or_else(|_| path.to_string());
    match std::fs::metadata(path) {
        Ok(m) => {
            let mtime = m
                .modified()
                .ok()
                .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
                .map(|d| d.as_nanos())
                .unwrap_or(0);
            format!("{canon}\t{}\t{mtime}", m.len())
        }
        Err(_) => format!("{canon}\tabsent"),
    }
}

/// The running executable's fingerprint: a rebuild changes size or mtime, so cached intermediates never
/// survive a code change.
pub fn exe_fingerprint() -> String {
    std::env::current_exe().map(|p| file_fingerprint(&p.display().to_string())).unwrap_or_else(|_| "unknown".into())
}

/// Every `RUSTLE_*` variable, sorted, one `k=v` per line, minus `exclude`.
pub fn env_fingerprint(exclude: &[&str]) -> String {
    let mut v: Vec<(String, String)> = std::env::vars()
        .filter(|(k, _)| k.starts_with("RUSTLE_") && !exclude.contains(&k.as_str()))
        .collect();
    v.sort();
    v.into_iter().map(|(k, val)| format!("env\t{k}={val}\n")).collect()
}

/// One keyed cache entry: `<root>/<kind>/<fnv(key)>/`.
pub struct Entry {
    pub dir: PathBuf,
    pub key: String,
    /// Payload files a complete entry of this kind must list in `DONE` (`reps`: reps.tsv + reps.fa, `paf`: out.paf).
    pub required: &'static [&'static str],
}

impl Entry {
    pub fn new(root: &Path, kind: &str, key: String) -> Entry {
        let dir = root.join(kind).join(format!("{:016x}", fnv1a64(key.as_bytes())));
        let required: &'static [&'static str] = match kind {
            "reps" => &["reps.tsv", "reps.fa"],
            "paf" => &["out.paf"],
            _ => &[],
        };
        Entry { dir, key, required }
    }
    /// A complete entry whose recorded key equals this one and whose files still have the sizes recorded at
    /// commit (`DONE` lists `name<TAB>bytes`), so a truncated file is a miss, not a silent partial replay.
    pub fn is_hit(&self) -> bool {
        let Ok(done) = std::fs::read_to_string(self.dir.join("DONE")) else { return false };
        // an empty or partial DONE (a crash between rename and write-back) is a miss, never a vacuous hit
        let listed: Vec<&str> = done.lines().filter_map(|l| l.split_once('\t').map(|(n, _)| n)).collect();
        if listed.is_empty() || !self.required.iter().all(|r| listed.contains(r)) {
            return false;
        }
        let sizes_ok = done.lines().all(|l| match l.split_once('\t') {
            Some((name, n)) => std::fs::metadata(self.dir.join(name)).map(|m| m.len().to_string() == n).unwrap_or(false),
            None => false,
        });
        sizes_ok && std::fs::read_to_string(self.dir.join("key.tsv")).map(|k| k == self.key).unwrap_or(false)
    }
    /// A fresh temporary directory next to the entry; [`Entry::commit`] renames it into place.
    pub fn staging(&self) -> Result<PathBuf> {
        static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let tmp = self.dir.with_extension(format!("tmp{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
        std::fs::write(tmp.join("key.tsv"), &self.key)?;
        Ok(tmp)
    }
    /// Durably publish a staging directory: every payload is fsynced, `DONE` (name + size of each file) is
    /// written and fsynced, the directory is renamed into place and the parent fsynced. If another run already
    /// published a complete entry for the same key, the staging copy is discarded instead.
    pub fn commit(&self, staging: &Path) -> Result<()> {
        if self.is_hit() {
            let _ = std::fs::remove_dir_all(staging);
            return Ok(());
        }
        let mut done = String::new();
        let mut names: Vec<String> = std::fs::read_dir(staging)?
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().to_string())
            .filter(|n| n != "DONE")
            .collect();
        names.sort();
        for n in names {
            std::fs::File::open(staging.join(&n))?.sync_all()?;
            done.push_str(&format!("{n}\t{}\n", std::fs::metadata(staging.join(&n))?.len()));
        }
        {
            let mut f = std::fs::File::create(staging.join("DONE"))?;
            f.write_all(done.as_bytes())?;
            f.sync_all()?;
        }
        let _ = std::fs::remove_dir_all(&self.dir); // an incomplete or stale entry under this name
        if let Some(p) = self.dir.parent() {
            std::fs::create_dir_all(p)?;
        }
        std::fs::rename(staging, &self.dir).with_context(|| format!("committing {}", self.dir.display()))?;
        if let Some(p) = self.dir.parent() {
            if let Ok(d) = std::fs::File::open(p) {
                let _ = d.sync_all();
            }
        }
        Ok(())
    }
}

/// Write the representatives: `reps.tsv` + `reps.fa`, exact and in index order.
pub fn write_reps(dir: &Path, reps: &[DenovoTranscript]) -> Result<()> {
    let mut t = std::io::BufWriter::new(std::fs::File::create(dir.join("reps.tsv"))?);
    writeln!(t, "idx\ttid\tchrom\tstart\tend\tn_reads\tstrand\tintrons\tdistinguishing_uniq\tcore_bp\tstub\ttes\tseq_len")?;
    let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join("reps.fa"))?);
    for (i, r) in reps.iter().enumerate() {
        let introns: Vec<String> = r.introns.iter().map(|(d, a)| format!("{d}-{a}")).collect();
        writeln!(
            t,
            "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.tid,
            r.chrom,
            r.start,
            r.end,
            r.n_reads,
            r.strand,
            if introns.is_empty() { "-".to_string() } else { introns.join(",") },
            r.distinguishing_uniq,
            r.core_bp,
            r.stub as u8,
            r.tes.map_or_else(|| "-".to_string(), |v| v.to_string()),
            r.seq.len()
        )?;
        writeln!(f, ">{i}")?;
        f.write_all(&r.seq)?;
        writeln!(f)?;
    }
    t.flush()?;
    f.flush()?;
    Ok(())
}

/// Read back what [`write_reps`] wrote.
pub fn read_reps(dir: &Path) -> Result<Vec<DenovoTranscript>> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    for (k, line) in std::io::BufReader::new(std::fs::File::open(dir.join("reps.fa"))?).split(b'\n').enumerate() {
        let line = line?;
        if k % 2 == 1 {
            seqs.push(line);
        }
    }
    let mut out = Vec::new();
    for (k, line) in std::io::BufReader::new(std::fs::File::open(dir.join("reps.tsv"))?).lines().enumerate() {
        let line = line?;
        if k == 0 {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        anyhow::ensure!(c.len() == 13, "reps.tsv: bad row {k}");
        let i: usize = c[0].parse()?;
        let introns: Vec<(u64, u64)> = if c[7] == "-" {
            Vec::new()
        } else {
            c[7].split(',')
                .map(|p| {
                    let (d, a) = p.split_once('-').context("intron")?;
                    Ok((d.parse()?, a.parse()?))
                })
                .collect::<Result<_>>()?
        };
        let seq = seqs.get(i).cloned().context("reps.fa shorter than reps.tsv")?;
        anyhow::ensure!(seq.len() == c[12].parse::<usize>()?, "reps.fa/reps.tsv length mismatch at {i}");
        out.push(DenovoTranscript {
            tid: c[1].to_string(),
            chrom: c[2].to_string(),
            start: c[3].parse()?,
            end: c[4].parse()?,
            n_reads: c[5].parse()?,
            strand: c[6].chars().next().context("strand")?,
            introns,
            seq,
            distinguishing_uniq: c[8].parse()?,
            core_bp: c[9].parse()?,
            stub: c[10] == "1",
            tes: if c[11] == "-" { None } else { Some(c[11].parse()?) },
        });
    }
    anyhow::ensure!(out.len() == seqs.len(), "reps.tsv/reps.fa row count mismatch");
    Ok(out)
}

/// `minimap2 --version`, once per process (part of every PAF key).
pub fn minimap2_version(minimap2: &str) -> String {
    static V: std::sync::OnceLock<std::sync::Mutex<std::collections::HashMap<String, String>>> =
        std::sync::OnceLock::new();
    let m = V.get_or_init(Default::default);
    if let Some(v) = m.lock().unwrap().get(minimap2) {
        return v.clone();
    }
    let v = std::process::Command::new(minimap2)
        .arg("--version")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|_| "unknown".into());
    m.lock().unwrap().insert(minimap2.to_string(), v.clone());
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(i: u64, introns: Vec<(u64, u64)>, seq: &[u8]) -> DenovoTranscript {
        DenovoTranscript {
            tid: format!("DN_c_{i}"),
            chrom: "c1".into(),
            start: 100 * i,
            end: 100 * i + 50,
            n_reads: i as u32 + 2,
            strand: if i % 2 == 0 { '+' } else { '-' },
            introns,
            seq: seq.to_vec(),
            distinguishing_uniq: i as usize,
            core_bp: 7 * i,
            stub: i % 3 == 0,
            tes: if i % 2 == 0 { Some(100 * i + 60) } else { None },
        }
    }

    #[test]
    fn reps_round_trip_exactly() {
        let dir = std::env::temp_dir().join(format!("rustle_run_cache_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let reps = vec![tx(0, vec![], b"ACGTacgtNN"), tx(1, vec![(110, 120), (130, 140)], b"GGGccc"), tx(3, vec![], b"")];
        write_reps(&dir, &reps).unwrap();
        let back = read_reps(&dir).unwrap();
        assert_eq!(back.len(), reps.len());
        for (a, b) in reps.iter().zip(&back) {
            assert_eq!(
                (&a.tid, &a.chrom, a.start, a.end, a.n_reads, a.strand, &a.introns, &a.seq, a.distinguishing_uniq, a.core_bp, a.stub, a.tes),
                (&b.tid, &b.chrom, b.start, b.end, b.n_reads, b.strand, &b.introns, &b.seq, b.distinguishing_uniq, b.core_bp, b.stub, b.tes)
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_hit_needs_done_and_the_identical_key() {
        let root = std::env::temp_dir().join(format!("rustle_run_cache_hit_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let e = Entry::new(&root, "reps", "key one\n".into());
        assert!(!e.is_hit());
        let st = e.staging().unwrap();
        assert!(!e.is_hit(), "staging is not a hit");
        write_reps(&st, &[]).unwrap(); // a reps entry needs its payload (reps.tsv + reps.fa) to be complete
        e.commit(&st).unwrap();
        assert!(e.is_hit());
        // same directory name forced, different key text: must miss
        let other = Entry { dir: e.dir.clone(), key: "key two\n".into(), required: e.required };
        assert!(!other.is_hit());
        // a file truncated after commit: must miss
        let e2 = Entry::new(&root, "paf", "k\n".into());
        let st2 = e2.staging().unwrap();
        std::fs::write(st2.join("out.paf"), b"line one\nline two\n").unwrap();
        e2.commit(&st2).unwrap();
        assert!(e2.is_hit());
        std::fs::write(e2.dir.join("out.paf"), b"line one\n").unwrap();
        assert!(!e2.is_hit(), "a truncated file must not replay");
        // an empty DONE (crash between rename and write-back) is a miss, not a vacuous hit
        std::fs::write(e2.dir.join("out.paf"), b"line one\nline two\n").unwrap();
        assert!(e2.is_hit());
        std::fs::write(e2.dir.join("DONE"), b"").unwrap();
        assert!(!e2.is_hit(), "an empty DONE must not be a hit");
        // a reps entry whose DONE lacks reps.fa is a miss
        let e3 = Entry::new(&root, "reps", "r\n".into());
        let st3 = e3.staging().unwrap();
        std::fs::write(st3.join("reps.tsv"), b"idx\n").unwrap();
        e3.commit(&st3).unwrap();
        assert!(!e3.is_hit(), "a reps entry without reps.fa must not be a hit");
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn fnv_incremental_matches_one_shot() {
        let mut f = Fnv::default();
        f.update(b"abc");
        f.update(b"def");
        assert_eq!(f.finish(), fnv1a64(b"abcdef"));
    }
}
