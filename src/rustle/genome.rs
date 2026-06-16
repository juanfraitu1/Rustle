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
    fn rt_switch_rejects_n_runs_and_edges() {
        // All-N matching flanks must NOT count as a repeat.
        let seq = b"NNNNNNNNCCCCCCCCNNNNNNNN";
        let g = GenomeIndex::from_seqs(&[("c1", seq)]);
        assert!(!g.is_rt_switch("c1", 8, 24, 4), "N-runs are not real repeats");
        // donor < repeat_len -> false (no room)
        assert!(!g.is_rt_switch("c1", 2, 24, 4));
    }
}
