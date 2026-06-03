//! Minimal FASTA loader and splice consensus checks for junction validation.

use crate::types::DetHashMap as HashMap;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};

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
        let mut entries: Vec<(String, usize, u64)> = Vec::new();
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
            if length == 0 {
                continue;
            }
            entries.push((f[0].to_string(), length, offset));
        }
        if entries.is_empty() {
            // Nothing to subset (no overlap) -> safest to load the whole genome.
            return Self::from_fasta(path);
        }
        let mut file = File::open(path).with_context(|| format!("failed to open FASTA: {}", path))?;
        let mut seqs: HashMap<String, Vec<u8>> = Default::default();
        for (name, length, offset) in entries {
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

    #[inline]
    fn base(&self, chrom: &str, pos0: u64) -> Option<u8> {
        let seq = self.seqs.get(chrom)?;
        seq.get(pos0 as usize).copied()
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

    /// Check reference-like splice consensus at intron boundaries.
    /// Junction coordinates are donor=left exon end, acceptor=right exon start (0-based half-open).
    /// The intron spans [donor, acceptor) — first intron base is at `donor`, last at `acceptor - 1`.
    pub fn is_consensus_splice(
        &self,
        chrom: &str,
        donor: u64,
        acceptor: u64,
        strand: Option<i8>,
    ) -> bool {
        // First two bases of intron: positions donor, donor+1 (0-based)
        let plus_left_1 = self.base(chrom, donor);
        let plus_left_2 = self.base(chrom, donor.saturating_add(1));
        // Last two bases of intron: positions acceptor-2, acceptor-1 (0-based)
        let plus_right_2 = acceptor.checked_sub(2).and_then(|p| self.base(chrom, p));
        let plus_right_1 = acceptor.checked_sub(1).and_then(|p| self.base(chrom, p));

        let plus_ok = matches!(plus_left_1, Some(b'G'))
            && matches!(plus_left_2, Some(b'T') | Some(b'C'))
            && matches!(plus_right_2, Some(b'A'))
            && matches!(plus_right_1, Some(b'G'));

        let minus_ok = matches!(plus_left_1, Some(b'C'))
            && matches!(plus_left_2, Some(b'T'))
            && matches!(plus_right_1, Some(b'C'))
            && matches!(plus_right_2, Some(b'A') | Some(b'G'));

        match strand {
            Some(1) => plus_ok,
            Some(-1) => minus_ok,
            _ => plus_ok || minus_ok,
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
}
