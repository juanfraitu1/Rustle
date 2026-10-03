//! On-disk cache of the catalog builder's expensive intermediates, for fast re-runs and for analysis.
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `RUSTLE_CACHE_DIR`, set by default by `tools/rustle_pipeline.sh`)
//!
//! Enabled by `RUSTLE_CACHE_DIR=<dir>` (the pipeline driver sets `PREFIX.cache`); unset = nothing is read or
//! written and every output is byte-identical to a build without this module.
//!
//! These objects are cached, each at a boundary where everything downstream reads only the cached object:
//!
//! * **`reps/<key>/`** — the collapsed locus representatives with their read statistics, i.e. the state of
//!   `detect_homology_catalog_genome_wide` after both BAM passes and the locus collapse (human chr16: ~260 of
//!   ~370 s). Files: `reps.tsv` (one row per representative, in representative-index order), `reps.fa` (the
//!   exact sequence bytes, one line each), `key.tsv` (the full key material), `DONE`.
//!   Under `gw_family_catalog --piecewise` the same kind holds one entry PER CONTIG (key line
//!   `rustle catalog reps v1 contig=<name>`) and one for their merge (`rustle catalog reps v1 piecewise-merge`),
//!   same files; see `denovo_pipeline::detect_homology_catalog_piecewise`. A contig split by `--piece-records` has
//!   one entry per read-free sub-range (`contig=<name>:<lo>-<hi>`) and its cut plan in **`plan/<key>/pieces.tsv`**.
//! * **`paf/<key>/`** — one all-vs-all minimap2 PAF per (query bytes, command line, minimap2 version):
//!   `out.paf`, `key.tsv`, `DONE`.
//!   The families stage (`mcl_families --from-gtf`, key line `rustle families paf v2`) keys its entry on a
//!   [`ContentHash`] of the loci FASTA taken WHILE the FASTA is written (no re-read), and replays a hit by
//!   HARD-LINKING `out.paf` to `<out>.loci.paf` ([`Entry::replay`]; a copy where a link is impossible) instead of
//!   copying it. Its entry is PINNED ([`Entry::pinned`]): `DONE` also records the payload's mtime (ns), inode, a
//!   sampled-content fingerprint and its full content hash, so an in-place write to the shared inode (through the
//!   linked product) is a miss, never a stale replay. The one theoretical stale replay left is the `.asbin` sidecar's
//!   edge: a same-size in-place rewrite that restores the mtime to the nanosecond and changes no sampled block;
//!   `RUSTLE_CACHE_VERIFY=1` re-hashes every pinned payload on a hit and closes it (audit mode). Every writer of a
//!   replayed product unlinks it first, never truncates it, and a payload is linked to at most one product (a
//!   second prefix sharing the cache directory gets a copy), so two products never share an inode.
//!   `o3_candidates::minimap2` (key line `rustle o3 minimap2 v1`) keeps each of its minimap2 calls here the same way, keyed on the
//!   content hashes of the target and of the query file.
//! * **`cand/<key>/`** — the result of the `o3_candidates` stage (spec `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`
//!   §5.8): `candidates.tsv` and `contigs.fa` are required (an empty `contigs.fa`, a run that flagged nothing, is a complete payload).
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
    "RUSTLE_CACHE_VERIFY",
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

/// A stable, word-at-a-time 128-bit content hash: two multiply-rotate lanes over the input's little-endian 8-byte
/// words (the last partial word zero-padded), then the byte length, then a murmur3 finaliser per lane. The value
/// depends only on the byte stream, never on how [`ContentHash::update`] calls split it, and not on the Rust release
/// or the platform. About 8x faster than byte-wise [`Fnv`]; used where a cache key must cover every byte of a large
/// file this process writes itself (hashed as it is written, so a hit never re-reads it).
#[derive(Clone, Debug)]
pub struct ContentHash {
    a: u64,
    b: u64,
    len: u64,
    tail: [u8; 8],
    ntail: usize,
}
impl Default for ContentHash {
    fn default() -> Self {
        ContentHash { a: 0x243f_6a88_85a3_08d3, b: 0x1319_8a2e_0370_7344, len: 0, tail: [0; 8], ntail: 0 }
    }
}
impl ContentHash {
    #[inline]
    fn word(&mut self, w: u64) {
        self.a = (self.a ^ w).wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(31);
        self.b = (self.b ^ w.rotate_left(23)).wrapping_mul(0xc2b2_ae3d_27d4_eb4f).rotate_left(27);
    }
    pub fn update(&mut self, mut bytes: &[u8]) {
        self.len += bytes.len() as u64;
        if self.ntail > 0 {
            let take = (8 - self.ntail).min(bytes.len());
            self.tail[self.ntail..self.ntail + take].copy_from_slice(&bytes[..take]);
            self.ntail += take;
            bytes = &bytes[take..];
            if self.ntail < 8 {
                return;
            }
            let w = u64::from_le_bytes(self.tail);
            self.word(w);
            self.ntail = 0;
        }
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            self.word(u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]));
        }
        let rem = chunks.remainder();
        self.tail[..rem.len()].copy_from_slice(rem);
        self.ntail = rem.len();
    }
    /// Bytes hashed so far.
    pub fn len(&self) -> u64 {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    /// 32 hex digits (the two lanes).
    pub fn hex(&self) -> String {
        fn fmix(mut k: u64) -> u64 {
            k ^= k >> 33;
            k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
            k ^= k >> 33;
            k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
            k ^ (k >> 33)
        }
        let mut s = self.clone();
        if s.ntail > 0 {
            let mut t = [0u8; 8];
            t[..s.ntail].copy_from_slice(&s.tail[..s.ntail]);
            s.word(u64::from_le_bytes(t));
        }
        s.word(s.len);
        format!("{:016x}{:016x}", fmix(s.a), fmix(s.b ^ s.a.rotate_left(32)))
    }
    /// The hash of a whole file, streamed in 4 MiB pieces.
    pub fn of_file(path: &Path) -> std::io::Result<ContentHash> {
        use std::io::Read;
        let mut f = std::fs::File::open(path)?;
        let mut h = ContentHash::default();
        let mut buf = vec![0u8; 4 << 20];
        loop {
            let n = f.read(&mut buf)?;
            if n == 0 {
                return Ok(h);
            }
            h.update(&buf[..n]);
        }
    }
}

/// A writer that hashes exactly the bytes it passes on (a [`ContentHash`] of the file as written), or passes them
/// on untouched when built with `enabled = false` (`hash` is then `None`).
pub struct HashingWriter<W: Write> {
    pub inner: W,
    pub hash: Option<ContentHash>,
}
impl<W: Write> HashingWriter<W> {
    pub fn new(inner: W, enabled: bool) -> Self {
        HashingWriter { inner, hash: enabled.then(ContentHash::default) }
    }
}
impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        if let Some(h) = self.hash.as_mut() {
            h.update(&buf[..n]);
        }
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// The pins of a file, `mtime_ns<TAB>inode<TAB>sample`: what an in-place write through a hard link changes (a link
/// does not). `sample` is an FNV-1a of 16 evenly spaced 4 KiB blocks (the last included) with their offsets, the
/// `.asbin` sidecar's sampling (`denovo_assemble::AsTsvIdentity`), so even a same-size rewrite with the mtime put
/// back is caught when it touches a sampled block.
fn file_pins(path: &Path) -> Option<String> {
    use std::io::{Read, Seek, SeekFrom};
    let m = std::fs::metadata(path).ok()?;
    #[cfg(unix)]
    let (t, ino) = {
        use std::os::unix::fs::MetadataExt;
        (m.mtime() as i128 * 1_000_000_000 + m.mtime_nsec() as i128, m.ino())
    };
    #[cfg(not(unix))]
    let (t, ino) = (m.modified().ok()?.duration_since(std::time::UNIX_EPOCH).ok()?.as_nanos() as i128, 0u64);
    let size = m.len();
    let mut f = std::fs::File::open(path).ok()?;
    let mut fnv = Fnv::default();
    let block = 4096u64;
    let mut buf: Vec<u8> = Vec::with_capacity(block as usize);
    for i in 0..16u64 {
        let at = if size <= block { 0 } else { (size - block) * i / 15 };
        f.seek(SeekFrom::Start(at)).ok()?;
        buf.clear();
        (&mut f).take(block).read_to_end(&mut buf).ok()?;
        fnv.update(&at.to_le_bytes());
        fnv.update(&buf);
    }
    Some(format!("{t}\t{ino}\t{:016x}", fnv.finish()))
}

/// `RUSTLE_CACHE_VERIFY=1`: re-hash every pinned payload on a hit (audit mode; see the module doc).
pub fn verify_mode() -> bool {
    std::env::var("RUSTLE_CACHE_VERIFY").map_or(false, |v| v == "1")
}

/// Link `src` to `dest` (a new name for the same inode) after unlinking `dest`, or copy when the two are on
/// different file systems, links are unsupported, or (`sole`) `src` already has another name besides its cache
/// entry. `dest` is never truncated in place, so an older inode that `dest` named (a cache payload it was linked to)
/// is left intact. Returns whether a link was made.
///
/// `sole` keeps a payload linked to at most ONE product: when several output prefixes share one cache directory,
/// the first replay links and every other prefix gets a copy, so an in-place write to one product can change the
/// cache entry (which its pins then turn into a miss) but never another prefix's product. The driver gives each
/// prefix its own `PREFIX.cache`, so a re-run there always links (it unlinks its own old product first).
pub fn link_or_copy(src: &Path, dest: &Path, sole: bool) -> std::io::Result<bool> {
    match std::fs::remove_file(dest) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e),
    }
    #[cfg(unix)]
    let unshared = {
        use std::os::unix::fs::MetadataExt;
        !sole || std::fs::metadata(src)?.nlink() == 1
    };
    #[cfg(not(unix))]
    let unshared = !sole;
    if unshared && std::fs::hard_link(src, dest).is_ok() {
        return Ok(true);
    }
    std::fs::copy(src, dest)?;
    Ok(false)
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
    /// Payload files a complete entry of this kind must list in `DONE` (`reps`: reps.tsv + reps.fa, `paf`: out.paf, `cand`: candidates.tsv + contigs.fa).
    pub required: &'static [&'static str],
    /// A PINNED entry's payloads may be hard-linked out ([`Entry::replay`], [`Entry::stage_link`]): its `DONE` lines
    /// are `name<TAB>bytes<TAB>mtime_ns<TAB>inode<TAB>sample<TAB>content_hash` and a hit needs all of them to match
    /// (the full content hash only under [`verify_mode`]), so a write through a linked product invalidates the entry.
    pub pin: bool,
}

impl Entry {
    pub fn new(root: &Path, kind: &str, key: String) -> Entry {
        let dir = root.join(kind).join(format!("{:016x}", fnv1a64(key.as_bytes())));
        let required: &'static [&'static str] = match kind {
            "reps" => &["reps.tsv", "reps.fa"],
            "paf" => &["out.paf"],
            "plan" => &["pieces.tsv"],
            "cand" => &["candidates.tsv", "contigs.fa"],
            _ => &[],
        };
        Entry { dir, key, required, pin: false }
    }
    /// This entry, pinned (see [`Entry::pin`]).
    pub fn pinned(mut self) -> Entry {
        self.pin = true;
        self
    }
    /// A complete entry whose recorded key equals this one and whose files still have the sizes recorded at
    /// commit (`DONE` lists `name<TAB>bytes`), so a truncated file is a miss, not a silent partial replay. A pinned
    /// entry's files must also keep the mtime and inode recorded at commit ([`Entry::pin`]).
    pub fn is_hit(&self) -> bool {
        self.is_hit_verify(self.pin && verify_mode())
    }
    /// [`Entry::is_hit`], with the payload re-hash of a pinned entry forced on or off (`verify`).
    pub fn is_hit_verify(&self, verify: bool) -> bool {
        let Ok(done) = std::fs::read_to_string(self.dir.join("DONE")) else { return false };
        // an empty or partial DONE (a crash between rename and write-back) is a miss, never a vacuous hit
        let listed: Vec<&str> = done.lines().filter_map(|l| l.split_once('\t').map(|(n, _)| n)).collect();
        if listed.is_empty() || !self.required.iter().all(|r| listed.contains(r)) {
            return false;
        }
        let sizes_ok = done.lines().all(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let path = self.dir.join(c[0]);
            let size_ok = c.len() >= 2
                && std::fs::metadata(&path).map(|m| m.len().to_string() == c[1]).unwrap_or(false);
            if !size_ok {
                return false;
            }
            if c.len() == 2 {
                return !self.pin; // a pinned entry must carry its pins
            }
            // pinned: `name bytes mtime_ns inode sample hash` (the pins of the file as committed)
            c.len() == 6
                && file_pins(&path).map_or(false, |p| p == c[2..5].join("\t"))
                && (!verify || ContentHash::of_file(&path).map_or(false, |h| h.hex() == c[5]))
        });
        sizes_ok && std::fs::read_to_string(self.dir.join("key.tsv")).map(|k| k == self.key).unwrap_or(false)
    }
    /// Replay payload `name` of a hit as `dest` without copying it: `dest` becomes a hard link to the cached file,
    /// after unlinking whatever `dest` was; a copy where a link is impossible or the payload is already linked to
    /// another product ([`link_or_copy`] `sole`). Only for a pinned entry, whose `DONE` pins turn any later in-place
    /// write through `dest` into a miss.
    pub fn replay(&self, name: &str, dest: &Path) -> std::io::Result<bool> {
        assert!(self.pin, "only a pinned entry's payload may be linked out");
        link_or_copy(&self.dir.join(name), dest, true)
    }
    /// Put the product `src` (just written) into a staging directory as payload `name` by hard link (a copy where
    /// impossible); see [`Entry::replay`].
    pub fn stage_link(&self, staging: &Path, name: &str, src: &Path) -> std::io::Result<bool> {
        assert!(self.pin, "only a pinned entry's payload may be linked in");
        link_or_copy(src, &staging.join(name), true)
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
            let path = staging.join(&n);
            std::fs::File::open(&path)?.sync_all()?;
            let len = std::fs::metadata(&path)?.len();
            if self.pin {
                let pins = file_pins(&path).context("pinning a cache payload")?;
                let h = ContentHash::of_file(&path)?.hex();
                done.push_str(&format!("{n}\t{len}\t{pins}\t{h}\n"));
            } else {
                done.push_str(&format!("{n}\t{len}\n"));
            }
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
        let other = Entry { dir: e.dir.clone(), key: "key two\n".into(), required: e.required, pin: false };
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

    /// The `cand` kind (the `o3_candidates` result) requires `candidates.tsv` and `contigs.fa`: a half-written entry is a miss, and an
    /// EMPTY `contigs.fa` (a run that flagged no candidate) is a complete payload, not a missing one.
    #[test]
    fn a_cand_entry_is_a_hit_only_with_both_its_candidates_table_and_its_contigs() {
        let dir = tempfile::tempdir().unwrap();
        let e = Entry::new(dir.path(), "cand", "rustle o3 candidates v1\n".into());
        assert_eq!(e.required, ["candidates.tsv", "contigs.fa"]);
        assert!(e.dir.starts_with(dir.path().join("cand")), "{}", e.dir.display());
        let st = e.staging().unwrap();
        std::fs::write(st.join("candidates.tsv"), b"family\n").unwrap();
        e.commit(&st).unwrap();
        assert!(!e.is_hit(), "a cand entry without contigs.fa must not be a hit");
        let st = e.staging().unwrap();
        std::fs::write(st.join("candidates.tsv"), b"family\n").unwrap();
        std::fs::write(st.join("contigs.fa"), b"").unwrap();
        e.commit(&st).unwrap();
        assert!(e.is_hit(), "both payloads present (one of them empty): a hit");
    }

    #[test]
    fn content_hash_is_the_byte_stream_whatever_the_split() {
        let data: Vec<u8> = (0..1000u32).map(|i| (i.wrapping_mul(2_654_435_761) >> 13) as u8).collect();
        let mut one = ContentHash::default();
        one.update(&data);
        for split in [&[1usize, 7, 8, 9, 100][..], &[3, 3, 3, 3], &[999], &[0, 0, 5, 0]] {
            let mut h = ContentHash::default();
            let mut at = 0;
            for &n in split {
                let n = n.min(data.len() - at);
                h.update(&data[at..at + n]);
                at += n;
            }
            h.update(&data[at..]);
            assert_eq!((h.hex(), h.len()), (one.hex(), 1000));
        }
        // zero padding of the last word is not ambiguous (the length is hashed), every byte counts
        let hx = |b: &[u8]| {
            let mut h = ContentHash::default();
            h.update(b);
            h.hex()
        };
        assert_ne!(hx(b"ab"), hx(b"ab\0"));
        assert_ne!(hx(b""), hx(b"\0"));
        let mut flipped = data.clone();
        for i in [0usize, 7, 8, 500, 999] {
            flipped[i] ^= 1;
            assert_ne!(hx(&flipped), one.hex(), "byte {i}");
            flipped[i] ^= 1;
        }
        // a HashingWriter hashes exactly what it writes; disabled, it only passes the bytes on
        let mut w = HashingWriter::new(Vec::new(), true);
        w.write_all(&data[..10]).unwrap();
        w.write_all(&data[10..]).unwrap();
        assert_eq!((w.inner.as_slice(), w.hash.unwrap().hex()), (&data[..], one.hex()));
        let mut off = HashingWriter::new(Vec::new(), false);
        off.write_all(&data).unwrap();
        assert!(off.hash.is_none() && off.inner == data);
    }

    /// A pinned entry's payload is replayed as a hard link; `dest` is unlinked, never truncated; a write through the
    /// link is a miss (mtime), a same-size rewrite with the mtime put back is a miss when it touches a sampled block,
    /// and the one edge left (outside every sampled block) is caught by the verify (full re-hash) mode.
    #[test]
    fn a_pinned_entry_replays_by_link_and_a_write_through_the_link_is_a_miss() {
        let root = std::env::temp_dir().join(format!("rustle_run_cache_pin_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let data: Vec<u8> = (0..(1u32 << 20)).map(|i| b"ACGT\t\n"[(i % 6) as usize]).collect();
        let product = root.join("run.loci.paf");
        std::fs::write(&product, &data).unwrap();
        let e = Entry::new(&root.join("cache"), "paf", "families k\n".into()).pinned();
        let st = e.staging().unwrap();
        assert!(e.stage_link(&st, "out.paf", &product).unwrap(), "same file system: a link, not a copy");
        e.commit(&st).unwrap();
        assert!(e.is_hit_verify(false) && e.is_hit_verify(true));
        let done = std::fs::read_to_string(e.dir.join("DONE")).unwrap();
        assert!(done.lines().all(|l| l.split('\t').count() == 6), "pinned DONE rows: {done}");
        // the same entry seen unpinned-style (a DONE without pins) must not satisfy a pinned lookup
        let plain: String = done.lines().map(|l| l.split('\t').take(2).collect::<Vec<_>>().join("\t") + "\n").collect();
        std::fs::write(e.dir.join("DONE"), &plain).unwrap();
        assert!(!e.is_hit_verify(false), "a pinned entry without pins is a miss");
        std::fs::write(e.dir.join("DONE"), &done).unwrap();
        assert!(e.is_hit_verify(false));
        #[cfg(unix)]
        let ino = |p: &Path| {
            use std::os::unix::fs::MetadataExt;
            std::fs::metadata(p).unwrap().ino()
        };
        #[cfg(unix)]
        assert_eq!(ino(&product), ino(&e.dir.join("out.paf")), "the product and the payload are one inode");
        // a SECOND prefix sharing the cache, whose old product names another inode: that inode is left intact, and
        // the payload (already linked to `product`) is copied, not linked, so the two products never alias
        let other = root.join("again.loci.paf");
        std::fs::write(&other, b"old product").unwrap();
        let keep = root.join("keep");
        std::fs::hard_link(&other, &keep).unwrap();
        assert!(!e.replay("out.paf", &other).unwrap(), "already linked to another product: a copy");
        assert_eq!(std::fs::read(&keep).unwrap(), b"old product");
        assert_eq!(std::fs::read(&other).unwrap(), data);
        #[cfg(unix)]
        assert_ne!(ino(&other), ino(&e.dir.join("out.paf")));
        // the same prefix re-run: its own product is unlinked first, so the payload is unshared again and linked
        let dest = product.clone();
        assert!(e.replay("out.paf", &dest).unwrap(), "a re-run of the linked prefix links again");
        assert_eq!(std::fs::read(&dest).unwrap(), data);
        #[cfg(unix)]
        assert_eq!(ino(&dest), ino(&e.dir.join("out.paf")));
        assert!(e.is_hit_verify(true), "a link changes neither mtime nor content");
        // an in-place same-size write through the link (a later run of an older binary, a shell redirect)
        let mtime = std::fs::metadata(&dest).unwrap().modified().unwrap();
        std::thread::sleep(std::time::Duration::from_millis(30));
        let mut w = data.clone();
        w[0] = b'X';
        std::fs::write(&dest, &w).unwrap();
        assert!(!e.is_hit_verify(false), "the mtime moved");
        // the same write with the mtime put back: caught by the sampled block at offset 0
        std::fs::File::options().write(true).open(&dest).unwrap().set_modified(mtime).unwrap();
        assert!(!e.is_hit_verify(false), "a sampled block changed");
        // the documented edge: a change outside every sampled block, mtime restored -> only the audit mode sees it
        let mut w2 = data.clone();
        w2[100_000] = b'X';
        std::fs::write(&dest, &w2).unwrap();
        std::fs::File::options().write(true).open(&dest).unwrap().set_modified(mtime).unwrap();
        assert!(e.is_hit_verify(false), "the documented edge (same size, mtime restored, no sampled block touched)");
        assert!(!e.is_hit_verify(true), "RUSTLE_CACHE_VERIFY=1 re-hashes the payload and rejects it");
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
