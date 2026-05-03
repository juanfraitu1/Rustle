//! Minimal FASTA loader and splice consensus checks for junction validation.

use crate::types::DetHashMap as HashMap;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};

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
