use crate::types::{DetHashMap as HashMap, DetHashSet};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::os::unix::fs::FileExt;

#[derive(Debug, Clone, Default)]
pub struct GenomeIndex {
    seqs: HashMap<String, Vec<u8>>,
    /// [`GenomeIndex::from_fasta_regions`]: per contig its length and its loaded windows `(start, bases)`, sorted and
    /// disjoint. A window-only contig answers [`GenomeIndex::fetch_sequence`] inside one window and is not listed by
    /// [`GenomeIndex::chroms`].
    windows: HashMap<String, (usize, Vec<(usize, Vec<u8>)>)>,
}

impl GenomeIndex {
    pub fn from_fasta(path: &str) -> Result<Self> {
        let f = File::open(path).with_context(|| format!("failed to open FASTA: {}", path))?;
        let reader = BufReader::new(f);
        let mut seqs: HashMap<String, Vec<u8>> = Default::default();
        let mut current: Option<String> = None;
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            if let Some(h) = line.strip_prefix('>') {
                let name = h.split_whitespace().next().unwrap_or("").to_string();
                if name.is_empty() {
                    current = None;
                    continue;
                }
                seqs.entry(name.clone()).or_default();
                current = Some(name);
            } else if let Some(chr) = current.as_ref() {
                let seq = seqs.get_mut(chr).expect("header inserted");
                seq.extend(
                    line.as_bytes()
                        .iter()
                        .filter(|c| !c.is_ascii_whitespace())
                        .map(|c| c.to_ascii_uppercase()),
                );
            }
        }
        Ok(Self {
            seqs,
            windows: Default::default(),
        })
    }

    /// A genome with no contigs loaded — for callers that provably fetch nothing (the catalog's cache-hit path
    /// when no re-admission option is on). `from_fasta_contigs` with an EMPTY set loads the whole genome instead.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Load ONLY the named contigs, seeking to each via the FASTA `.fai` index so
    /// the rest of the (possibly multi-GB) genome is never read. This is the fast
    /// path for region-scoped --vg runs. Falls back to a full `from_fasta` load if
    /// the `.fai` is missing, is malformed, or none of the wanted contigs are
    /// indexed — so the result is never worse than loading everything.
    pub fn from_fasta_contigs(path: &str, wanted: &DetHashSet<String>) -> Result<Self> {
        let fai = match std::fs::read_to_string(format!("{}.fai", path)) {
            Ok(s) => s,
            Err(_) => return Self::from_fasta(path),
        };
        // .fai columns: name \t length \t offset \t linebases \t linewidth
        let mut entries: Vec<(String, usize, u64, usize, usize)> = Vec::new();
        for line in fai.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            if !wanted.contains(f[0]) {
                continue;
            }
            let length: usize = match f[1].parse() {
                Ok(v) => v,
                Err(_) => return Self::from_fasta(path),
            };
            let offset: u64 = match f[2].parse() {
                Ok(v) => v,
                Err(_) => return Self::from_fasta(path),
            };
            let linebases: usize = f[3].parse().unwrap_or(0);
            let linewidth: usize = f[4].parse().unwrap_or(0);
            if length == 0 {
                continue;
            }
            entries.push((f[0].to_string(), length, offset, linebases, linewidth));
        }
        if entries.is_empty() {
            // Nothing to subset (no overlap) -> safest to load the whole genome.
            return Self::from_fasta(path);
        }
        let mut file =
            File::open(path).with_context(|| format!("failed to open FASTA: {}", path))?;
        // Contigs are independent: read them in parallel through positioned reads on one shared handle.
        let bulk: Vec<Option<Vec<u8>>> = {
            use rayon::prelude::*;
            let shared = &file;
            entries
                .par_iter()
                .map(|(_, length, offset, linebases, linewidth)| {
                    Self::read_contig_bulk(shared, *length, *offset, *linebases, *linewidth)
                })
                .collect::<Result<_>>()?
        };
        let mut seqs: HashMap<String, Vec<u8>> = Default::default();
        for ((name, length, offset, _, _), bulk) in entries.into_iter().zip(bulk) {
            if let Some(seq) = bulk {
                seqs.insert(name, seq);
                continue;
            }
            // Fallback: the byte-at-a-time reader. Reached when the `.fai` line geometry is absent,
            // non-uniform, or does not describe the bytes actually on disk.
            file.seek(SeekFrom::Start(offset))?;
            let mut reader = BufReader::new(&mut file);
            let mut seq: Vec<u8> = Vec::with_capacity(length);
            let mut lbuf = Vec::new();
            while seq.len() < length {
                lbuf.clear();
                let n = reader.read_until(b'\n', &mut lbuf)?;
                if n == 0 {
                    break;
                }
                for &c in &lbuf {
                    if !c.is_ascii_whitespace() {
                        seq.push(c.to_ascii_uppercase());
                        if seq.len() >= length {
                            break;
                        }
                    }
                }
            }
            seqs.insert(name, seq);
        }
        Ok(Self {
            seqs,
            windows: Default::default(),
        })
    }

    /// Read one contig as a SINGLE bulk read plus per-line `extend_from_slice`, using the `.fai` line
    /// geometry to know exactly where the terminators sit, then one `make_ascii_uppercase` over the result.
    ///
    /// Why: this loader was 22.8% of the 25-region control panel and is CPU-bound, not I/O-bound (~365
    /// Mbp/s), because the previous path called `read_until` per 60-base line and then tested and copied
    /// every byte individually. Measured 2.0x on real contigs (RABL2/GSTM/ACTB, two repeats each) with the
    /// produced bytes hash-identical to the byte-at-a-time path.
    ///
    /// Returns `Ok(None)` — meaning "caller must use the slow path" — whenever the geometry is not the
    /// uniform `linebases`/`linewidth` layout samtools writes, or whenever the bytes read do not match it.
    /// EVERY assumption is CHECKED rather than assumed: the terminator bytes must really be whitespace and
    /// the assembled sequence must be exactly `length` bases. This function must never return a sequence
    /// that differs from what the slow path would have produced — a wrong base here is silent, since
    /// `fetch_sequence` has no way to notice.
    /// Only the given regions (contig, 0-based start, end; merged per contig, clipped at the contig end), read in parallel
    /// through the `.fai`: the families stage needs its loci spans, a third of a genome or less, never whole contigs.
    /// [`GenomeIndex::fetch_sequence`] answers inside a region exactly as the whole-contig index does and None outside;
    /// contigs the `.fai` does not list are not loaded. Without a usable `.fai` (or line geometry) it is the regions'
    /// contigs' whole index ([`GenomeIndex::from_fasta_contigs`]).
    pub fn from_fasta_regions(path: &str, regions: &[(String, u64, u64)]) -> Result<Self> {
        let contigs: DetHashSet<String> = regions.iter().map(|r| r.0.clone()).collect();
        let Ok(fai) = std::fs::read_to_string(format!("{}.fai", path)) else {
            return Self::from_fasta_contigs(path, &contigs);
        };
        let mut geom: HashMap<String, (usize, u64, usize, usize)> = Default::default();
        for line in fai.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 || !contigs.contains(f[0]) {
                continue;
            }
            match (f[1].parse(), f[2].parse(), f[3].parse(), f[4].parse()) {
                (Ok(l), Ok(o), Ok(lb), Ok(lw)) => {
                    geom.insert(f[0].to_string(), (l, o, lb, lw));
                }
                _ => return Self::from_fasta_contigs(path, &contigs),
            }
        }
        // merge each contig's regions (overlapping or abutting) after clipping at its end
        let mut by: HashMap<String, Vec<(usize, usize)>> = Default::default();
        for (c, s, e) in regions {
            if let Some(&(len, ..)) = geom.get(c) {
                let (s, e) = ((*s as usize).min(len), (*e as usize).min(len));
                if s < e {
                    by.entry(c.clone()).or_default().push((s, e));
                }
            }
        }
        let mut jobs: Vec<(String, usize, usize)> = Vec::new();
        for (c, mut v) in by {
            v.sort_unstable();
            let mut cur: Option<(usize, usize)> = None;
            for (s, e) in v {
                cur = match cur {
                    Some((cs, ce)) if s <= ce => Some((cs, ce.max(e))),
                    Some(w) => {
                        jobs.push((c.clone(), w.0, w.1));
                        Some((s, e))
                    }
                    None => Some((s, e)),
                };
            }
            if let Some(w) = cur {
                jobs.push((c.clone(), w.0, w.1));
            }
        }
        let file = File::open(path).with_context(|| format!("failed to open FASTA: {}", path))?;
        let loaded: Vec<Option<Vec<u8>>> = {
            use rayon::prelude::*;
            let shared = &file;
            jobs.par_iter()
                .map(|(c, ws, we)| {
                    let (len, off, lb, lw) = geom[c];
                    Self::read_window(shared, len, off, lb, lw, *ws, *we)
                })
                .collect::<Result<_>>()?
        };
        if loaded.iter().any(|w| w.is_none()) {
            return Self::from_fasta_contigs(path, &contigs); // irregular line geometry: the careful whole-contig reader
        }
        let mut windows: HashMap<String, (usize, Vec<(usize, Vec<u8>)>)> = Default::default();
        for ((c, ws, _), w) in jobs.into_iter().zip(loaded) {
            let len = geom[&c].0;
            windows
                .entry(c)
                .or_insert_with(|| (len, Vec::new()))
                .1
                .push((ws, w.expect("checked above")));
        }
        for (_, wins) in windows.values_mut() {
            wins.sort_unstable_by_key(|w| w.0);
        }
        Ok(Self {
            seqs: Default::default(),
            windows,
        })
    }

    /// Bases `[ws, we)` of a contig with uniform `.fai` geometry: its whole lines are read and checked exactly as
    /// [`GenomeIndex::read_contig_bulk`] checks a whole contig, then trimmed; None when the geometry does not hold.
    fn read_window(
        file: &File,
        length: usize,
        offset: u64,
        linebases: usize,
        linewidth: usize,
        ws: usize,
        we: usize,
    ) -> Result<Option<Vec<u8>>> {
        if linebases == 0 || ws >= we || we > length {
            return Ok(None);
        }
        let (l0, l1) = (ws / linebases, (we - 1) / linebases + 1);
        Ok(Self::read_lines(file, length, offset, linebases, linewidth, l0, l1)?.map(|mut v| {
            v.truncate(we - l0 * linebases);
            v.drain(..ws - l0 * linebases);
            v
        }))
    }

    fn read_contig_bulk(
        file: &File,
        length: usize,
        offset: u64,
        linebases: usize,
        linewidth: usize,
    ) -> Result<Option<Vec<u8>>> {
        if linebases == 0 {
            return Ok(None);
        }
        let n_lines = length / linebases + usize::from(length % linebases > 0);
        let seq = Self::read_lines(file, length, offset, linebases, linewidth, 0, n_lines)?;
        Ok(seq.filter(|s| s.len() == length))
    }

    /// The bases of lines `[l0, l1)` of a contig with uniform geometry, uppercased; None when a line holds a terminator
    /// byte, a terminator is not whitespace, or the file is short (the last line of the contig may lack its terminator).
    fn read_lines(
        file: &File,
        length: usize,
        offset: u64,
        linebases: usize,
        linewidth: usize,
        l0: usize,
        l1: usize,
    ) -> Result<Option<Vec<u8>>> {
        // Uniform-geometry guard. `linewidth - linebases` is the terminator width ("\n" or "\r\n").
        if linebases == 0 || linewidth < linebases || linewidth - linebases > 2 {
            return Ok(None);
        }
        let term = linewidth - linebases;
        let n_full = length / linebases;
        let rem = length % linebases;
        let n_lines = n_full + usize::from(rem > 0);
        if l1 > n_lines || l0 > l1 {
            return Ok(None);
        }
        // Stream whole lines in blocks (positioned reads: the handle is shared across threads); the final line's
        // terminator may be absent (EOF), so a block may come back short.
        const LINES_PER_BLOCK: usize = 1 << 16;
        let mut buf = vec![0u8; LINES_PER_BLOCK.min((l1 - l0).max(1)) * linewidth];
        let mut seq: Vec<u8> = Vec::with_capacity((l1 - l0) * linebases);
        let (mut li, mut at) = (l0, offset + (l0 * linewidth) as u64);
        while li < l1 {
            let k = (l1 - li).min(LINES_PER_BLOCK);
            let want = k * linewidth;
            let mut got = 0usize;
            while got < want {
                let n = file.read_at(&mut buf[got..want], at + got as u64)?;
                if n == 0 {
                    break;
                }
                got += n;
            }
            let block = &buf[..got];
            let mut pos = 0usize;
            for j in 0..k {
                let line = li + j;
                let this = if line + 1 == n_lines && rem > 0 {
                    rem
                } else {
                    linebases
                };
                if pos + this > block.len() {
                    return Ok(None); // short file vs .fai — let the slow path decide
                }
                let chunk = &block[pos..pos + this];
                // The bases themselves must not contain a terminator; if they do the geometry is wrong.
                if chunk.iter().any(|c| c.is_ascii_whitespace()) {
                    return Ok(None);
                }
                seq.extend_from_slice(chunk);
                pos += this;
                // Verify the bytes we are about to SKIP really are line terminators.
                for t in 0..term {
                    match block.get(pos + t) {
                        Some(c) if c.is_ascii_whitespace() => {}
                        // Missing terminator is only acceptable at EOF on the last line.
                        None if line + 1 == n_lines => {}
                        _ => return Ok(None),
                    }
                }
                pos += term;
            }
            li += k;
            at += want as u64;
        }
        seq.make_ascii_uppercase();
        Ok(Some(seq))
    }

    /// Fetch a subsequence from the genome (0-based half-open coordinates).
    pub fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        if let Some(seq) = self.seqs.get(chrom) {
            let s = start as usize;
            let e = (end as usize).min(seq.len());
            if s >= e {
                return None;
            }
            return Some(seq[s..e].to_vec());
        }
        // a window-only contig: clip at the contig end as above, then the request must lie inside one window
        let (len, wins) = self.windows.get(chrom)?;
        let s = start as usize;
        let e = (end as usize).min(*len);
        if s >= e {
            return None;
        }
        let i = wins.partition_point(|(ws, _)| *ws <= s);
        let (ws, w) = wins.get(i.checked_sub(1)?)?;
        (e <= ws + w.len()).then(|| w[s - ws..e - ws].to_vec())
    }

    /// Iterate over (chrom_name, sequence_bytes) for every contig in the genome.
    /// Used by genome-wide scans (e.g. positional k-mer scan).
    pub fn chroms(&self) -> impl Iterator<Item = (&str, &[u8])> {
        self.seqs.iter().map(|(k, v)| (k.as_str(), v.as_slice()))
    }

    /// Length of a chromosome's sequence, or 0 if not present.
    pub fn chrom_len(&self, chrom: &str) -> u64 {
        self.seqs.get(chrom).map(|s| s.len() as u64).unwrap_or(0)
    }

    /// Canonical splice site (annotation-free, FASTA only). `donor` = left-exon-end,
    /// `acceptor` = right-exon-start (0-based half-open); the intron donor dinucleotide is
    /// `seq[donor..donor+2]`, the acceptor dinucleotide is `seq[acceptor-2..acceptor]`.
    /// Accepts the three canonical classes (and their reverse complements on '−'):
    ///   + strand: GT-AG (U2 major), GC-AG (minor), AT-AC (U12 minor)
    ///   − strand: CT-AC,            CT-GC,          GT-AT  (revcomp of the above)
    /// strand '.' accepts either orientation. Including GC-AG/AT-AC matters: GT-AG-only
    /// over-rejects real minor-class introns (the read-coherence realness gate was
    /// dropping genuine GC-AG isoforms).
    pub fn is_canonical_junction(
        &self,
        chrom: &str,
        donor: u64,
        acceptor: u64,
        strand: char,
    ) -> bool {
        if acceptor < donor + 4 {
            return false;
        }
        let d = match self.fetch_sequence(chrom, donor, donor + 2) {
            Some(s) => s,
            None => return false,
        };
        let a = match self.fetch_sequence(chrom, acceptor - 2, acceptor) {
            Some(s) => s,
            None => return false,
        };
        if d.len() < 2 || a.len() < 2 {
            return false;
        }
        let up = |b: &[u8]| -> [u8; 2] { [b[0].to_ascii_uppercase(), b[1].to_ascii_uppercase()] };
        let (d, a) = (up(&d), up(&a));
        let plus = matches!(
            (d, a),
            ([b'G', b'T'], [b'A', b'G'])  // GT-AG  (U2 major)
                | ([b'G', b'C'], [b'A', b'G'])  // GC-AG  (minor)
                | ([b'A', b'T'], [b'A', b'C']) // AT-AC  (U12 minor)
        );
        let minus = matches!(
            (d, a),
            ([b'C', b'T'], [b'A', b'C'])  // revcomp GT-AG
                | ([b'C', b'T'], [b'G', b'C'])  // revcomp GC-AG
                | ([b'G', b'T'], [b'A', b'T']) // revcomp AT-AC
        );
        match strand {
            '+' => plus,
            '-' => minus,
            _ => plus || minus,
        }
    }

    /// RT template-switch signature: a direct repeat flanking the junction (the
    /// `repeat_len` bp ending at the donor equals the `repeat_len` bp ending at the
    /// acceptor). Heuristic; N-runs do not count. Annotation-free (FASTA only).
    pub fn is_rt_switch(&self, chrom: &str, donor: u64, acceptor: u64, repeat_len: u64) -> bool {
        if donor < repeat_len || acceptor < repeat_len {
            return false;
        }
        let up = match self.fetch_sequence(chrom, donor - repeat_len, donor) {
            Some(s) => s,
            None => return false,
        };
        let dn = match self.fetch_sequence(chrom, acceptor - repeat_len, acceptor) {
            Some(s) => s,
            None => return false,
        };
        !up.is_empty()
            && up.eq_ignore_ascii_case(&dn)
            && !up.iter().any(|&b| b == b'N' || b == b'n')
    }

    /// Microhomology (direct-repeat) at a recombination breakpoint bracket `(left, right)` — the
    /// template-switch signature applied to an arbitrary switch point, not only a splice junction.
    /// Returns true if ANY repeat length in `[k_min, k_max]` shows a direct repeat (the `is_rt_switch`
    /// test) whose matched window is NOT low-complexity. A confirmed mosaic that carries this signature
    /// is more likely an RT/template-switch artifact than a biological gene conversion. Widening beyond
    /// the old fixed 8 bp (the documented `is_rt_switch` limitation) raises sensitivity to hotspots.
    ///
    /// LOW-COMPLEXITY GUARD: a homopolymer / dinucleotide-repeat window (e.g. `AAAAAA`, `ATATAT`)
    /// trivially matches a direct repeat almost everywhere, so without this filter a true gene
    /// conversion near a simple repeat would be wrongly demoted to an artifact (the error direction
    /// SUPPRESSES real conversions). We require the matched window to carry ≥ 3 distinct bases, which
    /// keeps informative repeats (`CGTACGTA`) and rejects the uninformative low-complexity ones.
    pub fn breakpoint_microhomology(
        &self,
        chrom: &str,
        left: u64,
        right: u64,
        k_min: u64,
        k_max: u64,
    ) -> bool {
        (k_min..=k_max).any(|k| {
            if !self.is_rt_switch(chrom, left, right, k) || left < k {
                return false;
            }
            // is_rt_switch already proved left-window == right-window; check the left window's complexity.
            match self.fetch_sequence(chrom, left - k, left) {
                Some(w) => !is_low_complexity_window(&w),
                None => false,
            }
        })
    }
}

/// A window is low-complexity (uninformative for a direct-repeat call) if it carries fewer than 3
/// distinct bases — i.e. a homopolymer (1) or a dinucleotide repeat (2). Such windows match a direct
/// repeat almost everywhere, so they must not trigger a microhomology call.
fn is_low_complexity_window(w: &[u8]) -> bool {
    let mut seen = [false; 4];
    let mut distinct = 0usize;
    for &b in w {
        let idx = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => continue,
        };
        if !seen[idx] {
            seen[idx] = true;
            distinct += 1;
        }
    }
    distinct < 3
}

// ===========================================================================
// Case-carrying indexed FASTA (pysam.FastaFile.fetch equivalent)
// ===========================================================================

#[derive(Clone, Copy, Debug)]
struct FaiRecord {
    length: u64,
    offset: u64,
    linebases: u64,
    linewidth: u64,
}

/// Minimal `.fai`-indexed FASTA reader returning bytes VERBATIM (soft-mask lowercase
/// preserved), matching `pysam.FastaFile.fetch(chrom, start, end)` on 0-based
/// half-open coordinates. Needed because `GenomeIndex` upper-cases and so
/// cannot feed a soft-mask computation.
pub struct IndexedFasta {
    file: File,
    index: std::collections::HashMap<String, FaiRecord>,
}

impl IndexedFasta {
    /// Open `path` and read its `path.fai` sidecar
    /// (`name \t length \t offset \t linebases \t linewidth`).
    pub fn open(path: &str) -> std::io::Result<Self> {
        let fai_text = std::fs::read_to_string(format!("{}.fai", path))?;
        let mut index = std::collections::HashMap::new();
        for line in fai_text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            let rec = FaiRecord {
                length: f[1].parse().unwrap_or(0),
                offset: f[2].parse().unwrap_or(0),
                linebases: f[3].parse().unwrap_or(0),
                linewidth: f[4].parse().unwrap_or(0),
            };
            if rec.linebases == 0 || rec.linewidth == 0 {
                continue;
            }
            index.insert(f[0].to_string(), rec);
        }
        let file = File::open(path)?;
        Ok(Self { file, index })
    }

    /// Fetch `[start, end)` (0-based half-open) verbatim, or `None` if the contig is
    /// unknown. An empty/degenerate interval yields `Some(vec![])` (pysam returns "").
    /// Bytes are read straight from the file at their `.fai`-computed offsets, skipping
    /// the per-line newline, so lowercase soft-mask is preserved exactly.
    pub fn fetch(&self, chrom: &str, start: i64, end: i64) -> Option<Vec<u8>> {
        let rec = self.index.get(chrom)?;
        let mut s = start.max(0) as u64;
        let mut e = end.max(0) as u64;
        if e > rec.length {
            e = rec.length;
        }
        if s > rec.length {
            s = rec.length;
        }
        if s >= e {
            return Some(Vec::new());
        }
        let mut out = Vec::with_capacity((e - s) as usize);
        let (lb, lw) = (rec.linebases, rec.linewidth);
        let mut p = s;
        while p < e {
            let line = p / lb;
            let col = p % lb;
            let byte_off = rec.offset + line * lw + col;
            let take = std::cmp::min(lb - col, e - p); // bases left on this line
            let mut buf = vec![0u8; take as usize];
            self.file.read_exact_at(&mut buf, byte_off).ok()?;
            out.extend_from_slice(&buf);
            p += take;
        }
        Some(out)
    }
}

#[cfg(test)]
impl GenomeIndex {
    /// Test-only constructor from an in-memory (chrom -> sequence) map. Sequences
    /// are uppercased to match `from_fasta`.
    pub(crate) fn from_seqs(pairs: &[(&str, &[u8])]) -> Self {
        let mut seqs: HashMap<String, Vec<u8>> = Default::default();
        for (name, seq) in pairs {
            seqs.insert(
                name.to_string(),
                seq.iter().map(|c| c.to_ascii_uppercase()).collect(),
            );
        }
        Self {
            seqs,
            windows: Default::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_fasta_contigs_loads_only_wanted_via_fai() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa");
        // c1: 2 lines (10 + 4 bases); c2: 1 line (10 bases).
        std::fs::write(&fa, ">c1\nACGTACGTAC\nACGT\n>c2\nTTTTGGGGAA\n").unwrap();
        // .fai cols: name length offset linebases linewidth
        //  ">c1\n"=4 -> c1 data @4, 14 bases; c1 data = 11+5 = 16 bytes;
        //  ">c2\n" @20 -> c2 data @24, 10 bases.
        std::fs::write(
            format!("{}.fai", fa.display()),
            "c1\t14\t4\t10\t11\nc2\t10\t24\t10\t11\n",
        )
        .unwrap();
        let mut wanted = DetHashSet::default();
        wanted.insert("c1".to_string());
        let g = GenomeIndex::from_fasta_contigs(fa.to_str().unwrap(), &wanted).unwrap();
        assert_eq!(
            g.fetch_sequence("c1", 0, 14),
            Some(b"ACGTACGTACACGT".to_vec()),
            "c1 must load exactly via the .fai offset"
        );
        assert!(
            g.fetch_sequence("c2", 0, 10).is_none(),
            "c2 must NOT be loaded (not in wanted set)"
        );
    }

    /// Every contig the `.fai` lists loads to exactly the whole-file parse (uppercased, terminators stripped), across
    /// line widths, CRLF, a final line without terminator and contigs longer than one read block.
    #[test]
    fn from_fasta_contigs_equals_the_whole_file_parse_for_every_contig() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa");
        let mut text = String::new();
        let mut fai = String::new();
        let mut x: u64 = 7;
        let specs: [(&str, usize, usize, &str); 6] =
            [("a", 1, 60, "\n"), ("b", 59, 60, "\n"), ("c", 500_007, 7, "\n"), ("d", 1_000_003, 61, "\n"), ("e", 33, 10, "\r\n"), ("f", 25, 10, "\n")];
        for (i, (name, len, width, term)) in specs.iter().enumerate() {
            text.push_str(&format!(">{name} desc\n"));
            let offset = text.len();
            let seq: String = (0..*len)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGTacgtN"[(x >> 33) as usize % 9] as char
                })
                .collect();
            for (k, chunk) in seq.as_bytes().chunks(*width).enumerate() {
                text.push_str(std::str::from_utf8(chunk).unwrap());
                let last = (k + 1) * width >= *len;
                if !(last && i + 1 == specs.len()) {
                    text.push_str(term); // the very last line of the file has no terminator
                }
            }
            fai.push_str(&format!("{name}\t{len}\t{offset}\t{width}\t{}\n", width + term.len()));
        }
        std::fs::write(&fa, &text).unwrap();
        std::fs::write(format!("{}.fai", fa.display()), &fai).unwrap();
        let wanted: DetHashSet<String> = specs.iter().map(|s| s.0.to_string()).collect();
        let sub = GenomeIndex::from_fasta_contigs(fa.to_str().unwrap(), &wanted).unwrap();
        let whole = GenomeIndex::from_fasta(fa.to_str().unwrap()).unwrap();
        for (name, len, _, _) in specs {
            assert_eq!(sub.fetch_sequence(name, 0, len as u64), whole.fetch_sequence(name, 0, len as u64), "{name}");
            assert_eq!(sub.fetch_sequence(name, 0, len as u64).unwrap().len(), len, "{name}");
        }
    }

    /// A region index ([`GenomeIndex::from_fasta_regions`]) answers every fetch inside its (merged) regions exactly as the
    /// whole-contig index does, clips at the contig end the same way, and answers None outside its regions.
    #[test]
    fn from_fasta_regions_fetches_like_the_whole_contig_index_inside_its_regions() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa");
        let mut text = String::new();
        let mut fai = String::new();
        let mut x: u64 = 11;
        let specs: [(&str, usize, usize, &str); 4] =
            [("a", 500_007, 7, "\n"), ("b", 1_003, 60, "\r\n"), ("c", 61, 61, "\n"), ("d", 2_500, 50, "\n")];
        for (i, (name, len, width, term)) in specs.iter().enumerate() {
            text.push_str(&format!(">{name}\n"));
            let offset = text.len();
            let seq: String = (0..*len)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGTacgtN"[(x >> 33) as usize % 9] as char
                })
                .collect();
            for (k, chunk) in seq.as_bytes().chunks(*width).enumerate() {
                text.push_str(std::str::from_utf8(chunk).unwrap());
                if !((k + 1) * width >= *len && i + 1 == specs.len()) {
                    text.push_str(term);
                }
            }
            fai.push_str(&format!("{name}\t{len}\t{offset}\t{width}\t{}\n", width + term.len()));
        }
        std::fs::write(&fa, &text).unwrap();
        std::fs::write(format!("{}.fai", fa.display()), &fai).unwrap();
        let regions: Vec<(String, u64, u64)> = vec![
            ("a".into(), 0, 13), ("a".into(), 10, 700), ("a".into(), 250_000, 260_003), ("a".into(), 499_990, 500_007),
            ("b".into(), 59, 61), ("b".into(), 900, 2_000), ("c".into(), 0, 61), ("d".into(), 1_234, 1_235),
            ("zz".into(), 0, 10),
        ];
        let reg = GenomeIndex::from_fasta_regions(fa.to_str().unwrap(), &regions).unwrap();
        let all: DetHashSet<String> = specs.iter().map(|s| s.0.to_string()).collect();
        let whole = GenomeIndex::from_fasta_contigs(fa.to_str().unwrap(), &all).unwrap();
        // every sub-interval of every merged region, including the clipped one at b's end
        let merged: [(&str, u64, u64); 6] =
            [("a", 0, 700), ("a", 250_000, 260_003), ("a", 499_990, 500_007), ("b", 59, 61), ("b", 900, 1_003), ("c", 0, 61)];
        for (c, s0, e0) in merged {
            for s in [s0, s0 + 1, (s0 + e0) / 2, e0 - 1] {
                for e in [s + 1, (s + e0).div_ceil(2), e0] {
                    if s < e {
                        assert_eq!(reg.fetch_sequence(c, s, e), whole.fetch_sequence(c, s, e), "{c}:{s}-{e}");
                    }
                }
            }
        }
        assert_eq!(reg.fetch_sequence("b", 950, 5_000), whole.fetch_sequence("b", 950, 5_000), "clipped at the contig end");
        assert_eq!(reg.fetch_sequence("d", 1_234, 1_235), whole.fetch_sequence("d", 1_234, 1_235));
        assert_eq!(reg.fetch_sequence("a", 690, 710), None, "crosses the region end");
        assert_eq!(reg.fetch_sequence("a", 1_000, 1_010), None, "outside every region");
        assert_eq!(reg.fetch_sequence("zz", 0, 5), None, "a contig the .fai does not list");
    }

    #[test]
    fn from_fasta_contigs_falls_back_to_full_without_fai() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa");
        std::fs::write(&fa, ">c1\nACGT\n>c2\nTTTT\n").unwrap(); // no .fai written
        let mut wanted = DetHashSet::default();
        wanted.insert("c1".to_string());
        let g = GenomeIndex::from_fasta_contigs(fa.to_str().unwrap(), &wanted).unwrap();
        // fallback = full load -> both contigs present
        assert_eq!(g.fetch_sequence("c1", 0, 4), Some(b"ACGT".to_vec()));
        assert_eq!(g.fetch_sequence("c2", 0, 4), Some(b"TTTT".to_vec()));
    }

    // is_canonical_junction: intron [donor,acceptor) is GT..AG on '+' (canonical),
    // CT..AC on '-' (the reverse complement). Build a sequence whose intron begins
    // GT and ends AG, and another that begins GG (non-canonical).
    #[test]
    fn canonical_junction_gt_ag_plus() {
        // positions:        0123456789...
        //   exon  [0,5)  = "AAAAA"
        //   intron[5,15) = "GT......AG"  (GT at 5,6; AG at 13,14)
        //   exon  [15,20)= "CCCCC"
        let seq = b"AAAAAGTTTTTTTAGCCCCC";
        let g = GenomeIndex::from_seqs(&[("c1", seq)]);
        // donor=5 (exon end), acceptor=15 (next exon start); intron = [5,15)
        assert!(
            g.is_canonical_junction("c1", 5, 15, '+'),
            "GT..AG on '+' is canonical"
        );
        // strand '.' accepts either orientation
        assert!(g.is_canonical_junction("c1", 5, 15, '.'));
        // wrong strand: GT..AG is NOT CT..AC, so '-' rejects it
        assert!(!g.is_canonical_junction("c1", 5, 15, '-'));
    }

    #[test]
    fn canonical_junction_ct_ac_minus() {
        // intron[5,15) = "CT......AC"
        let seq = b"AAAAACTTTTTTTACCCCCC";
        let g = GenomeIndex::from_seqs(&[("c1", seq)]);
        assert!(
            g.is_canonical_junction("c1", 5, 15, '-'),
            "CT..AC on '-' is canonical"
        );
        assert!(g.is_canonical_junction("c1", 5, 15, '.'));
        assert!(!g.is_canonical_junction("c1", 5, 15, '+'));
    }

    // Build a 20bp seq whose intron is [5,15): exon[0,5) + donor[5,7) + 6bp filler[7,13)
    // + acceptor[13,15) + exon[15,20). Constructed (not hand-counted) so donor/acceptor
    // dinucleotides land exactly where is_canonical_junction reads them.
    fn splice_seq(donor: &[u8; 2], acceptor: &[u8; 2]) -> Vec<u8> {
        let mut v = Vec::new();
        v.extend_from_slice(b"AAAAA"); // exon  [0,5)
        v.extend_from_slice(donor); //    intron[5,7)
        v.extend_from_slice(b"GGGGGG"); // filler[7,13)
        v.extend_from_slice(acceptor); // intron[13,15)
        v.extend_from_slice(b"CCCCC"); // exon  [15,20)
        v
    }

    #[test]
    fn canonical_junction_gc_ag_and_at_ac_plus() {
        // GC-AG minor class
        let gc = splice_seq(b"GC", b"AG");
        let g = GenomeIndex::from_seqs(&[("c1", &gc[..])]);
        assert!(
            g.is_canonical_junction("c1", 5, 15, '+'),
            "GC..AG on '+' is canonical (minor)"
        );
        assert!(g.is_canonical_junction("c1", 5, 15, '.'));
        assert!(
            !g.is_canonical_junction("c1", 5, 15, '-'),
            "GC..AG is not a '-' class"
        );
        // AT-AC U12 class
        let at = splice_seq(b"AT", b"AC");
        let g2 = GenomeIndex::from_seqs(&[("c1", &at[..])]);
        assert!(
            g2.is_canonical_junction("c1", 5, 15, '+'),
            "AT..AC on '+' is canonical (U12)"
        );
        assert!(!g2.is_canonical_junction("c1", 5, 15, '-'));
    }

    #[test]
    fn canonical_junction_minor_revcomp_minus() {
        // revcomp(GC-AG) = CT-GC
        let ctgc = splice_seq(b"CT", b"GC");
        let g = GenomeIndex::from_seqs(&[("c1", &ctgc[..])]);
        assert!(
            g.is_canonical_junction("c1", 5, 15, '-'),
            "CT..GC on '-' is canonical (revcomp GC-AG)"
        );
        assert!(!g.is_canonical_junction("c1", 5, 15, '+'));
        // revcomp(AT-AC) = GT-AT
        let gtat = splice_seq(b"GT", b"AT");
        let g2 = GenomeIndex::from_seqs(&[("c1", &gtat[..])]);
        assert!(
            g2.is_canonical_junction("c1", 5, 15, '-'),
            "GT..AT on '-' is canonical (revcomp AT-AC)"
        );
        // GT..AT must NOT pass on '+' (would need AG acceptor)
        assert!(!g2.is_canonical_junction("c1", 5, 15, '+'));
    }

    #[test]
    fn canonical_junction_noncanonical_rejected() {
        // intron[5,15) = "GG......AG" -> donor GG is non-canonical for both strands
        let seq = b"AAAAAGGTTTTTTAGCCCCC";
        let g = GenomeIndex::from_seqs(&[("c1", seq)]);
        assert!(!g.is_canonical_junction("c1", 5, 15, '+'));
        assert!(!g.is_canonical_junction("c1", 5, 15, '-'));
        assert!(!g.is_canonical_junction("c1", 5, 15, '.'));
    }

    #[test]
    fn canonical_junction_too_short_or_missing() {
        let seq = b"AAAAAGTAGCCCCC";
        let g = GenomeIndex::from_seqs(&[("c1", seq)]);
        // acceptor < donor+4 -> false (degenerate)
        assert!(!g.is_canonical_junction("c1", 5, 7, '+'));
        // missing chrom -> false
        assert!(!g.is_canonical_junction("nope", 5, 15, '+'));
    }

    // is_rt_switch: the `repeat_len` bp ending at the donor equals the `repeat_len`
    // bp ending at the acceptor (direct-repeat template-switch signature).
    #[test]
    fn rt_switch_direct_repeat_true() {
        // Make the 4 bp ending at donor (positions 4..8) equal the 4 bp ending at
        // acceptor (positions 16..20). Put "ACGT" at both [4,8) and [16,20).
        //          0123 4567 89... 16..20
        // index:   AAAA ACGT XXXXXXXX ACGT
        let mut seq = vec![b'A'; 24];
        seq[4..8].copy_from_slice(b"ACGT");
        seq[16..20].copy_from_slice(b"ACGT");
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        // donor=8 -> up = seq[4..8] = ACGT ; acceptor=20 -> dn = seq[16..20] = ACGT
        assert!(
            g.is_rt_switch("c1", 8, 20, 4),
            "matching flanks => RT-switch"
        );
    }

    #[test]
    fn rt_switch_no_repeat_false() {
        // Flanks differ.
        let mut seq = vec![b'A'; 24];
        seq[4..8].copy_from_slice(b"ACGT");
        seq[16..20].copy_from_slice(b"TGCA");
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        assert!(
            !g.is_rt_switch("c1", 8, 20, 4),
            "differing flanks => not RT-switch"
        );
    }

    #[test]
    fn breakpoint_microhomology_scans_k_range() {
        // An EXACTLY-4 bp direct repeat ([4,8)==[16,20)) bounded by DIFFERING flank bases (seq[3]!=seq[15])
        // so it does not extend: found when the k-range includes 4, missed when the range starts above 4.
        let mut seq = vec![b'A'; 24];
        seq[3] = b'G';
        seq[15] = b'T'; // distinct flanks just outside the repeat -> repeat length is exactly 4
        seq[4..8].copy_from_slice(b"CGTC");
        seq[16..20].copy_from_slice(b"CGTC");
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        assert!(
            g.breakpoint_microhomology("c1", 8, 20, 4, 6),
            "k-range covering 4 finds the repeat"
        );
        assert!(
            !g.breakpoint_microhomology("c1", 8, 20, 5, 6),
            "k-range above the 4bp repeat misses it"
        );
    }

    #[test]
    fn breakpoint_microhomology_rejects_low_complexity_window() {
        // Two poly-A loci: the windows ending at 100 and 200 are both "AAAA..." -> is_rt_switch would
        // call a direct repeat, but the low-complexity guard rejects it (homopolymer => 1 distinct base),
        // so a real conversion near a poly-A tract is NOT wrongly demoted to an RT-switch artifact.
        let seq = vec![b'A'; 300];
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        assert!(
            g.is_rt_switch("c1", 100, 200, 8),
            "poly-A trivially matches as a direct repeat"
        );
        assert!(
            !g.breakpoint_microhomology("c1", 100, 200, 6, 12),
            "low-complexity (homopolymer) window must NOT count as microhomology"
        );
    }

    /// The bulk contig reader must produce EXACTLY the bytes the byte-at-a-time reader produced, on the
    /// geometries it accepts (soft-masked bases, CRLF, a final line with no terminator, and a last line
    /// shorter than `linebases`) AND must fall back rather than guess when the geometry is not uniform.
    /// A wrong base here is silent — `fetch_sequence` cannot notice — so this is checked, not assumed.
    #[test]
    fn bulk_contig_reader_matches_the_byte_at_a_time_reader() {
        // (label, file bytes, fai line, expected sequence)
        let cases: [(&str, &[u8], &str, &[u8]); 5] = [
            (
                "LF, ragged last line",
                b">c1\nACGTACGTAC\nACGT\n",
                "c1\t14\t4\t10\t11\n",
                b"ACGTACGTACACGT",
            ),
            (
                "soft-masked bases must be uppercased",
                b">c1\nacgtACGTac\nnNnN\n",
                "c1\t14\t4\t10\t11\n",
                b"ACGTACGTACNNNN",
            ),
            (
                "no terminator on the final line",
                b">c1\nACGTACGTAC\nACGT",
                "c1\t14\t4\t10\t11\n",
                b"ACGTACGTACACGT",
            ),
            (
                "CRLF line terminators",
                b">c1\r\nACGTACGTAC\r\nACGT\r\n",
                "c1\t14\t5\t10\t12\n",
                b"ACGTACGTACACGT",
            ),
            (
                // Non-uniform geometry: the .fai claims 10 bases/line but line 1 holds 6. The bulk path
                // must DECLINE and let the slow path produce the (identical) answer.
                "non-uniform lines -> fallback",
                b">c1\nACGTAC\nGTACACGT\n",
                "c1\t14\t4\t10\t11\n",
                b"ACGTACGTACACGT",
            ),
        ];
        for (label, bytes, fai, expect) in cases {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("g.fa");
            std::fs::write(&fa, bytes).unwrap();
            std::fs::write(format!("{}.fai", fa.display()), fai).unwrap();
            let mut want = DetHashSet::default();
            want.insert("c1".to_string());
            let g = GenomeIndex::from_fasta_contigs(fa.to_str().unwrap(), &want).unwrap();
            assert_eq!(
                g.fetch_sequence("c1", 0, expect.len() as u64).unwrap(),
                expect.to_vec(),
                "bulk contig reader disagreed with the expected bytes: {label}"
            );
            assert_eq!(g.chrom_len("c1"), expect.len() as u64, "{label}");
        }
    }

    #[test]
    fn rt_switch_rejects_n_runs_and_edges() {
        // All-N matching flanks must NOT count as a repeat.
        let seq = b"NNNNNNNNCCCCCCCCNNNNNNNN";
        let g = GenomeIndex::from_seqs(&[("c1", seq)]);
        assert!(
            !g.is_rt_switch("c1", 8, 24, 4),
            "N-runs are not real repeats"
        );
        // donor < repeat_len -> false (no room)
        assert!(!g.is_rt_switch("c1", 2, 24, 4));
    }

    /// Case-carrying FASTA reader sanity on a tiny synthetic multi-line FASTA + `.fai`
    /// (validates the fai offset math + soft-mask/lowercase preservation without the
    /// 3.4G genome). Also exercises clamping + empty-interval behaviour.
    #[test]
    fn indexed_fasta_case_carrying() {
        let dir = std::env::temp_dir();
        let base = dir.join(format!("rustle_faidx_{}.fa", std::process::id()));
        let fa_path = base.to_str().unwrap().to_string();
        // 26 bases across 3 lines of linebases=10 (linewidth=11 with '\n'); mixed case.
        let seq = "ACGTacgtNN"; // line 1 (10)
        let seq2 = "GGGGcccctt"; // line 2 (10)
        let seq3 = "AAAAaa"; //     line 3 (6)
        std::fs::write(&fa_path, format!(">c1\n{seq}\n{seq2}\n{seq3}\n")).unwrap();
        // .fai: name length offset linebases linewidth ; offset = len(">c1\n") = 4
        std::fs::write(format!("{fa_path}.fai"), "c1\t26\t4\t10\t11\n").unwrap();

        let fa = IndexedFasta::open(&fa_path).unwrap();
        let whole: Vec<u8> = format!("{seq}{seq2}{seq3}").into_bytes();
        // full sequence, verbatim case
        assert_eq!(fa.fetch("c1", 0, 26).unwrap(), whole);
        // an interior slice straddling two lines, preserving lowercase
        assert_eq!(fa.fetch("c1", 4, 14).unwrap(), b"acgtNNGGGG".to_vec());
        // clamped end
        assert_eq!(fa.fetch("c1", 20, 999).unwrap(), b"AAAAaa".to_vec());
        // empty / unknown
        assert_eq!(fa.fetch("c1", 5, 5).unwrap(), Vec::<u8>::new());
        assert_eq!(fa.fetch("nope", 0, 3), None);

        let _ = std::fs::remove_file(&fa_path);
        let _ = std::fs::remove_file(format!("{fa_path}.fai"));
    }

    /// If the real genome is present, prove the case-carrying reader matches pysam by
    /// re-fetching 12 loci's exon sequences and comparing to the pysam bytes recorded in
    /// `vg_family/testdata/indexed_fasta_pysam_fixture.json` (trimmed from the retired
    /// `repeat_catalog_fixture.json`). Skipped (passes) when the 3.4G genome is unavailable (e.g. CI).
    #[test]
    fn indexed_fasta_matches_pysam_on_real_genome() {
        let genome = "/home/juanfra/winloci_scratch/GGO.fasta"; // test-only; not a shipped default
        if !std::path::Path::new(&format!("{genome}.fai")).exists() {
            eprintln!("real genome absent -> skipping pysam cross-check");
            return;
        }
        let fx: serde_json::Value =
            serde_json::from_str(include_str!("testdata/indexed_fasta_pysam_fixture.json"))
                .expect("parse indexed_fasta fixture json");
        let fa = match IndexedFasta::open(genome) {
            Ok(f) => f,
            Err(_) => return,
        };
        let mut n = 0usize;
        for locus in fx["loci"].as_array().unwrap() {
            let chrom = locus["chrom"].as_str().unwrap();
            let want = locus["seq"].as_str().unwrap().as_bytes();
            let mut got = Vec::new();
            for p in locus["exons"].as_array().unwrap() {
                let a = p.as_array().unwrap();
                got.extend(
                    fa.fetch(chrom, a[0].as_i64().unwrap(), a[1].as_i64().unwrap())
                        .unwrap(),
                );
            }
            assert_eq!(got, want, "pysam vs IndexedFasta MISMATCH at {chrom}");
            n += 1;
        }
        assert!(n >= 12, "fixture holds 12 loci, read {n}");
    }
}
