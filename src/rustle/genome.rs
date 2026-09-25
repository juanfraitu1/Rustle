//! Minimal FASTA loaders and splice consensus checks for junction validation.
//!
//! - `GenomeIndex` loads whole contigs (uppercased) into memory.
//! - `IndexedFasta` reads `[start, end)` straight from the file through its `.fai`, case preserved
//!   (soft-mask kept) — moved here from the removed `vg_family::repeat_catalog` (2026-09-24).

use crate::types::DetHashMap as HashMap;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::os::unix::fs::FileExt;

#[derive(Debug, Clone, Default)]
pub struct GenomeIndex {
    seqs: HashMap<String, Vec<u8>>,
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
        Ok(Self { seqs })
    }

    /// A genome with no contigs loaded — for callers that provably fetch nothing (the catalog's cache-hit path
    /// when no re-admission option is on). `from_fasta_contigs` with an EMPTY set loads the whole genome instead.
    pub fn empty() -> Self {
        Self { seqs: Default::default() }
    }

    /// Load ONLY the named contigs, seeking to each via the FASTA `.fai` index so
    /// the rest of the (possibly multi-GB) genome is never read. This is the fast
    /// path for region-scoped --vg runs. Falls back to a full `from_fasta` load if
    /// the `.fai` is missing, is malformed, or none of the wanted contigs are
    /// indexed — so the result is never worse than loading everything.
    pub fn from_fasta_contigs(
        path: &str,
        wanted: &std::collections::HashSet<String>,
    ) -> Result<Self> {
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
            let length: usize = match f[1].parse() { Ok(v) => v, Err(_) => return Self::from_fasta(path) };
            let offset: u64 = match f[2].parse() { Ok(v) => v, Err(_) => return Self::from_fasta(path) };
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
        let mut file = File::open(path).with_context(|| format!("failed to open FASTA: {}", path))?;
        let mut seqs: HashMap<String, Vec<u8>> = Default::default();
        for (name, length, offset, linebases, linewidth) in entries {
            if let Some(seq) =
                Self::read_contig_bulk(&mut file, length, offset, linebases, linewidth)?
            {
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
        Ok(Self { seqs })
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
    fn read_contig_bulk(
        file: &mut File,
        length: usize,
        offset: u64,
        linebases: usize,
        linewidth: usize,
    ) -> Result<Option<Vec<u8>>> {
        use std::io::Read;
        // Uniform-geometry guard. `linewidth - linebases` is the terminator width ("\n" or "\r\n").
        if linebases == 0 || linewidth < linebases || linewidth - linebases > 2 {
            return Ok(None);
        }
        let term = linewidth - linebases;
        let n_full = length / linebases;
        let rem = length % linebases;
        let n_lines = n_full + usize::from(rem > 0);
        // The final line's terminator may be absent (EOF), so ask for it but tolerate a short read.
        let want = length + n_lines * term;
        file.seek(SeekFrom::Start(offset))?;
        let mut buf = vec![0u8; want];
        let mut got = 0usize;
        while got < want {
            let n = file.read(&mut buf[got..])?;
            if n == 0 {
                break;
            }
            got += n;
        }
        buf.truncate(got);

        let mut seq: Vec<u8> = Vec::with_capacity(length);
        let mut pos = 0usize;
        for li in 0..n_lines {
            let this = if li + 1 == n_lines && rem > 0 { rem } else { linebases };
            if pos + this > buf.len() {
                return Ok(None); // short file vs .fai — let the slow path decide
            }
            let chunk = &buf[pos..pos + this];
            // The bases themselves must not contain a terminator; if they do the geometry is wrong.
            if chunk.iter().any(|c| c.is_ascii_whitespace()) {
                return Ok(None);
            }
            seq.extend_from_slice(chunk);
            pos += this;
            // Verify the bytes we are about to SKIP really are line terminators.
            for k in 0..term {
                match buf.get(pos + k) {
                    Some(c) if c.is_ascii_whitespace() => {}
                    // Missing terminator is only acceptable at EOF on the last line.
                    None if li + 1 == n_lines => {}
                    _ => return Ok(None),
                }
            }
            pos += term;
        }
        if seq.len() != length {
            return Ok(None);
        }
        seq.make_ascii_uppercase();
        Ok(Some(seq))
    }

    /// Fetch a subsequence from the genome (0-based half-open coordinates).
    pub fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        let seq = self.seqs.get(chrom)?;
        let s = start as usize;
        let e = (end as usize).min(seq.len());
        if s >= e {
            return None;
        }
        Some(seq[s..e].to_vec())
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
    pub fn is_canonical_junction(&self, chrom: &str, donor: u64, acceptor: u64, strand: char) -> bool {
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
    pub fn breakpoint_microhomology(&self, chrom: &str, left: u64, right: u64, k_min: u64, k_max: u64) -> bool {
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
        Self { seqs }
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
        let mut wanted = std::collections::HashSet::new();
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

    #[test]
    fn from_fasta_contigs_falls_back_to_full_without_fai() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa");
        std::fs::write(&fa, ">c1\nACGT\n>c2\nTTTT\n").unwrap(); // no .fai written
        let mut wanted = std::collections::HashSet::new();
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
        assert!(g.is_canonical_junction("c1", 5, 15, '+'), "GT..AG on '+' is canonical");
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
        assert!(g.is_canonical_junction("c1", 5, 15, '-'), "CT..AC on '-' is canonical");
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
        assert!(g.is_canonical_junction("c1", 5, 15, '+'), "GC..AG on '+' is canonical (minor)");
        assert!(g.is_canonical_junction("c1", 5, 15, '.'));
        assert!(!g.is_canonical_junction("c1", 5, 15, '-'), "GC..AG is not a '-' class");
        // AT-AC U12 class
        let at = splice_seq(b"AT", b"AC");
        let g2 = GenomeIndex::from_seqs(&[("c1", &at[..])]);
        assert!(g2.is_canonical_junction("c1", 5, 15, '+'), "AT..AC on '+' is canonical (U12)");
        assert!(!g2.is_canonical_junction("c1", 5, 15, '-'));
    }

    #[test]
    fn canonical_junction_minor_revcomp_minus() {
        // revcomp(GC-AG) = CT-GC
        let ctgc = splice_seq(b"CT", b"GC");
        let g = GenomeIndex::from_seqs(&[("c1", &ctgc[..])]);
        assert!(g.is_canonical_junction("c1", 5, 15, '-'), "CT..GC on '-' is canonical (revcomp GC-AG)");
        assert!(!g.is_canonical_junction("c1", 5, 15, '+'));
        // revcomp(AT-AC) = GT-AT
        let gtat = splice_seq(b"GT", b"AT");
        let g2 = GenomeIndex::from_seqs(&[("c1", &gtat[..])]);
        assert!(g2.is_canonical_junction("c1", 5, 15, '-'), "GT..AT on '-' is canonical (revcomp AT-AC)");
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
        assert!(g.is_rt_switch("c1", 8, 20, 4), "matching flanks => RT-switch");
    }

    #[test]
    fn rt_switch_no_repeat_false() {
        // Flanks differ.
        let mut seq = vec![b'A'; 24];
        seq[4..8].copy_from_slice(b"ACGT");
        seq[16..20].copy_from_slice(b"TGCA");
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        assert!(!g.is_rt_switch("c1", 8, 20, 4), "differing flanks => not RT-switch");
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
        assert!(g.breakpoint_microhomology("c1", 8, 20, 4, 6), "k-range covering 4 finds the repeat");
        assert!(!g.breakpoint_microhomology("c1", 8, 20, 5, 6), "k-range above the 4bp repeat misses it");
    }

    #[test]
    fn breakpoint_microhomology_rejects_low_complexity_window() {
        // Two poly-A loci: the windows ending at 100 and 200 are both "AAAA..." -> is_rt_switch would
        // call a direct repeat, but the low-complexity guard rejects it (homopolymer => 1 distinct base),
        // so a real conversion near a poly-A tract is NOT wrongly demoted to an RT-switch artifact.
        let seq = vec![b'A'; 300];
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        assert!(g.is_rt_switch("c1", 100, 200, 8), "poly-A trivially matches as a direct repeat");
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
            let mut want = std::collections::HashSet::new();
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
        assert!(!g.is_rt_switch("c1", 8, 24, 4), "N-runs are not real repeats");
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
            serde_json::from_str(include_str!("vg_family/testdata/indexed_fasta_pysam_fixture.json"))
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
                got.extend(fa.fetch(chrom, a[0].as_i64().unwrap(), a[1].as_i64().unwrap()).unwrap());
            }
            assert_eq!(got, want, "pysam vs IndexedFasta MISMATCH at {chrom}");
            n += 1;
        }
        assert!(n >= 12, "fixture holds 12 loci, read {n}");
    }
}
