//! Opt-in PSV-FASTA emission (`RUSTLE_VG_PSV_FASTA`).
//!
//! StringTie is sequence-blind: it emits coordinates, never the bases that distinguish paralog
//! copies. To *prove* Rustle's VG mode actually found copy-distinguishing variants (PSVs =
//! paralog-sequence-variants), this module emits, alongside the GTF:
//!   * `<out>.transcripts.fa` — the spliced nucleotide sequence of every emitted transcript
//!     (one record per GTF transcript, reverse-complemented for `-` strand), header carrying the
//!     family/copy id and PSV count.
//!   * `<out>.psv.tsv` — one row per copy-distinguishing site of each family transcript:
//!     `transcript_id  chrom  genomic_pos(1-based)  transcript_pos(1-based)  this_base
//!      sibling_bases  family_id  copy_id`.
//!
//! PSVs only exist inside multi-copy paralog families, so non-family transcripts get `n_psv=0` and
//! no `.psv.tsv` rows. Default-off → no files, GTF byte-identical.
//!
//! Caveat: copy diffs reuse `vg::diff_gf_exons`, which is column-wise within a syntenic exon pair, so
//! sites inside an indel-divergent block can be shifted/spurious; SNPs in well-anchored regions are
//! the high-confidence subset.

use std::collections::HashMap;
use std::io::{self, Write};

use crate::genome::GenomeIndex;
use crate::path_extract::Transcript;

/// One copy-distinguishing site at an absolute genomic position. Sites are joined to emitted
/// transcripts by GENOMIC OVERLAP (ref_pos ∈ transcript exons), NOT by `(family_id, copy_id)`:
/// that pair is only a handle to a parent BUNDLE, which assembly fans out into multiple
/// transcripts at different loci, so the copy_id on a transcript does not pin the locus its PSVs
/// live at. The `family_id`/`copy_id` here are the PSV's own provenance (the family-graph copy whose
/// exon carries the variant), reported for traceability.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PsvSite {
    /// Absolute genome-forward coordinate (0-based) of the variant.
    pub ref_pos: u64,
    /// The base this copy carries at `ref_pos`.
    pub this_base: u8,
    /// The distinct bases the sibling copies carry at `ref_pos` (sorted, deduped).
    pub sibling_bases: Vec<u8>,
    /// Provenance: the family this PSV's copy belongs to.
    pub family_id: usize,
    /// Provenance: the family-graph copy id (position in the family's bundle_indices) carrying it.
    pub copy_id: usize,
}

/// Build the spliced (mature) transcript sequence: concatenate the genome over the transcript's
/// exons in genomic-ascending order, then reverse-complement for a `-`-strand transcript so the
/// emitted sequence is in transcription order. Exons that fail to fetch contribute nothing.
pub fn spliced_sequence(
    genome: &GenomeIndex,
    chrom: &str,
    strand: char,
    exons: &[(u64, u64)],
) -> Vec<u8> {
    let mut sorted: Vec<(u64, u64)> = exons.to_vec();
    sorted.sort_unstable();
    let mut seq: Vec<u8> = Vec::new();
    for &(s, e) in &sorted {
        if let Some(chunk) = genome.fetch_sequence(chrom, s, e) {
            seq.extend_from_slice(&chunk);
        }
    }
    if strand == '-' {
        seq = crate::vg::reverse_complement(&seq);
    }
    seq
}

/// Map an absolute genomic position to a 0-based offset into the EMITTED (strand-oriented) spliced
/// sequence. `exons_sorted` must be genomic-ascending half-open `[s, e)`. Returns `None` if `gpos`
/// is intronic or outside the transcript (such a PSV is not present in the spliced sequence).
pub fn genomic_to_transcript_pos(
    exons_sorted: &[(u64, u64)],
    strand: char,
    gpos: u64,
) -> Option<u64> {
    let mut acc: u64 = 0;
    let mut total: u64 = 0;
    let mut fwd: Option<u64> = None;
    for &(s, e) in exons_sorted {
        let len = e.saturating_sub(s);
        if gpos >= s && gpos < e {
            fwd = Some(acc + (gpos - s));
        }
        acc += len;
        total += len;
    }
    let fwd = fwd?;
    if strand == '-' {
        // emitted sequence is the reverse-complement: index from the 3' (genomic-high) end.
        Some(total.saturating_sub(1).saturating_sub(fwd))
    } else {
        Some(fwd)
    }
}

/// One paralog copy within a family, as clustered from the EMITTED transcripts: a maximal set of
/// same-(chrom,strand) transcripts whose spans overlap. `rep_exons` are the genomic exons of the
/// most-complete member (most exons), used as the copy's representative geometry.
struct EmittedCopy {
    chrom: String,
    strand: char,
    rep_exons: Vec<(u64, u64)>,
}

/// Group emitted transcripts into per-family copies by genomic locus. Within each `vg_family_id`,
/// transcripts on the same chrom+strand whose spans overlap are one copy (isoforms); disjoint spans
/// are distinct copies. Non-family transcripts (no `vg_family_id`) are ignored. Returns
/// `family_id -> Vec<EmittedCopy>` (each copy's `rep_exons` = the member with the most exons).
fn cluster_emitted_copies(transcripts: &[Transcript]) -> HashMap<usize, Vec<EmittedCopy>> {
    let mut by_family: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, tx) in transcripts.iter().enumerate() {
        if let Some(f) = tx.vg_family_id {
            by_family.entry(f).or_default().push(i);
        }
    }
    let span = |tx: &Transcript| -> (u64, u64) {
        let lo = tx.exons.iter().map(|e| e.0).min().unwrap_or(0);
        let hi = tx.exons.iter().map(|e| e.1).max().unwrap_or(0);
        (lo, hi)
    };
    let mut out: HashMap<usize, Vec<EmittedCopy>> = HashMap::new();
    for (fam, mut idxs) in by_family {
        idxs.sort_by(|&a, &b| {
            let (ta, tb) = (&transcripts[a], &transcripts[b]);
            ta.chrom
                .cmp(&tb.chrom)
                .then(ta.strand.cmp(&tb.strand))
                .then(span(ta).0.cmp(&span(tb).0))
        });
        let mut copies: Vec<(String, char, u64, u64, usize)> = Vec::new(); // chrom,strand,lo,hi,rep_idx
        for &i in &idxs {
            let tx = &transcripts[i];
            let (lo, hi) = span(tx);
            match copies.last_mut() {
                Some((c, s, _clo, chi, rep))
                    if *c == tx.chrom && *s == tx.strand && lo < *chi =>
                {
                    // overlaps current copy → same copy; extend; keep most-complete representative.
                    if hi > *chi {
                        *chi = hi;
                    }
                    if tx.exons.len() > transcripts[*rep].exons.len() {
                        *rep = i;
                    }
                }
                _ => copies.push((tx.chrom.clone(), tx.strand, lo, hi, i)),
            }
        }
        let emitted: Vec<EmittedCopy> = copies
            .into_iter()
            .map(|(chrom, strand, _, _, rep)| {
                let mut rep_exons = transcripts[rep].exons.clone();
                rep_exons.sort_unstable();
                EmittedCopy { chrom, strand, rep_exons }
            })
            .collect();
        out.insert(fam, emitted);
    }
    out
}

/// Compute PSV sites DIRECTLY from the emitted transcripts, indexed by chromosome. For each family
/// (grouped by `vg_family_id`) the transcripts are clustered into copies by locus; same-strand copy
/// pairs are diffed via `vg::diagnostic_sites_between` over their representative exon sequences
/// (fetched genome-forward). Because the sites are computed from the emitted geometry, every PSV
/// ref_pos lies in a copy's representative exons — so the locus-overlap join in `emit_psv_fasta`
/// finds them (unlike bundle-derived sites, which sit at homology-shadow loci no transcript covers).
/// Per-chrom lists are sorted by ref_pos; sibling bases for the same (ref_pos, family, copy) are
/// merged. Opposite-strand (inverted) copy pairs are skipped (their forward sequences aren't
/// homologous).
pub fn build_psv_map(
    transcripts: &[Transcript],
    genome: &GenomeIndex,
) -> HashMap<String, Vec<PsvSite>> {
    let mut out: HashMap<String, Vec<PsvSite>> = HashMap::new();
    let fam_copies = cluster_emitted_copies(transcripts);
    // Deterministic family order.
    let mut fams: Vec<usize> = fam_copies.keys().copied().collect();
    fams.sort_unstable();
    let exon_seqs = |c: &EmittedCopy| -> Vec<(u64, Vec<u8>)> {
        c.rep_exons
            .iter()
            .map(|&(s, e)| (s, genome.fetch_sequence(&c.chrom, s, e).unwrap_or_default()))
            .collect()
    };
    for fam in fams {
        let copies = &fam_copies[&fam];
        if copies.len() < 2 {
            continue;
        }
        // Fetch each copy's representative exon sequences ONCE (not O(copies) times in the pair loop).
        let seqs: Vec<Vec<(u64, Vec<u8>)>> = copies.iter().map(&exon_seqs).collect();
        for a in 0..copies.len() {
            for b in (a + 1)..copies.len() {
                if copies[a].strand != copies[b].strand {
                    continue; // inverted pair: forward sequences not homologous
                }
                let sites = crate::vg::diagnostic_sites_between(&seqs[a], &seqs[b]);
                for (pa, ba, pb, bb) in sites {
                    out.entry(copies[a].chrom.clone()).or_default().push(PsvSite {
                        ref_pos: pa,
                        this_base: ba,
                        sibling_bases: vec![bb],
                        family_id: fam,
                        copy_id: a,
                    });
                    out.entry(copies[b].chrom.clone()).or_default().push(PsvSite {
                        ref_pos: pb,
                        this_base: bb,
                        sibling_bases: vec![ba],
                        family_id: fam,
                        copy_id: b,
                    });
                }
            }
        }
    }
    // Sort + merge sibling bases for the same (ref_pos, family, copy).
    for sites in out.values_mut() {
        sites.sort_by(|a, b| {
            a.ref_pos
                .cmp(&b.ref_pos)
                .then(a.family_id.cmp(&b.family_id))
                .then(a.copy_id.cmp(&b.copy_id))
                .then(a.this_base.cmp(&b.this_base))
        });
        let mut merged: Vec<PsvSite> = Vec::with_capacity(sites.len());
        for s in sites.drain(..) {
            match merged.last_mut() {
                Some(prev)
                    if prev.ref_pos == s.ref_pos
                        && prev.family_id == s.family_id
                        && prev.copy_id == s.copy_id =>
                {
                    for b in s.sibling_bases {
                        if !prev.sibling_bases.contains(&b) {
                            prev.sibling_bases.push(b);
                        }
                    }
                }
                _ => merged.push(s),
            }
        }
        for s in &mut merged {
            s.sibling_bases.sort_unstable();
            s.sibling_bases.dedup();
        }
        *sites = merged;
    }
    out
}

/// The transcript id Rustle writes in the GTF for transcript `i` (matches `gtf::write_gtf`):
/// the explicit `transcript_id` if present, else `{label}.{gene_no}.{tx_no}`.
fn transcript_id_for(tx: &Transcript, gno: usize, tno: usize, label: &str) -> String {
    if let (Some(_), Some(t)) = (&tx.gene_id, &tx.transcript_id) {
        t.clone()
    } else {
        format!("{}.{}.{}", label, gno.max(1), tno.max(1))
    }
}

/// Complement a single base (ACGT, case-preserving); other bytes pass through.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
        x => x,
    }
}

fn write_wrapped(w: &mut impl Write, seq: &[u8], width: usize) -> io::Result<()> {
    let mut i = 0;
    while i < seq.len() {
        let end = (i + width).min(seq.len());
        w.write_all(&seq[i..end])?;
        w.write_all(b"\n")?;
        i = end;
    }
    if seq.is_empty() {
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// Emit `<out_base>.transcripts.fa` (every transcript's spliced sequence) and `<out_base>.psv.tsv`
/// (one row per PSV whose absolute genomic position falls inside a transcript's exons). PSVs are
/// joined to transcripts by LOCUS OVERLAP (not by copy_id), so the result is correct even when a
/// paralog copy spans several bundles or a bundle fans out into transcripts at multiple loci.
/// Returns `(FASTA records, PSV rows)`.
pub fn emit_psv_fasta(
    transcripts: &[Transcript],
    genome: &GenomeIndex,
    label: &str,
    out_base: &std::path::Path,
) -> io::Result<(usize, usize)> {
    let psv_by_chrom = build_psv_map(transcripts, genome);
    let gene_tx_no = crate::gtf::assign_gene_tx_numbers(transcripts);
    let fa_path = out_base.with_extension("transcripts.fa");
    let tsv_path = out_base.with_extension("psv.tsv");
    let mut fa = io::BufWriter::new(std::fs::File::create(&fa_path)?);
    let mut tsv = io::BufWriter::new(std::fs::File::create(&tsv_path)?);
    writeln!(
        tsv,
        "transcript_id\tchrom\tgenomic_pos\ttranscript_pos\tthis_base\tsibling_bases\tfamily_id\tcopy_id"
    )?;

    let mut n_records = 0usize;
    let mut n_rows = 0usize;
    for (i, tx) in transcripts.iter().enumerate() {
        let (gno, tno) = gene_tx_no[i];
        let tid = transcript_id_for(tx, gno, tno, label);
        // Clamp exon ends to the contig length so spliced_sequence (which fetches clamped bytes) and
        // genomic_to_transcript_pos (which sums exon lengths) share ONE length model — otherwise an
        // exon running past the contig end (mismatched/truncated FASTA) would desync transcript_pos.
        let chrom_len = genome.chrom_len(&tx.chrom);
        let mut exons_sorted: Vec<(u64, u64)> = tx
            .exons
            .iter()
            .map(|&(s, e)| (s, e.min(chrom_len)))
            .filter(|&(s, e)| s < e)
            .collect();
        exons_sorted.sort_unstable();
        let seq = spliced_sequence(genome, &tx.chrom, tx.strand, &exons_sorted);

        // Locus-overlap join: PSVs on this chrom whose ref_pos lands inside an exon of THIS
        // transcript (intronic / out-of-range sites map to None and are skipped).
        let mut matched: Vec<(&PsvSite, u64)> = Vec::new();
        if let Some(sites) = psv_by_chrom.get(&tx.chrom) {
            for site in sites {
                if let Some(tpos) =
                    genomic_to_transcript_pos(&exons_sorted, tx.strand, site.ref_pos)
                {
                    matched.push((site, tpos));
                }
            }
        }
        let n_psv = matched.len();

        // FASTA record. family_id/copy_id is the transcript's own (bundle-level) vg tag, a label;
        // n_psv is the locus-accurate count of PSVs actually present in this emitted sequence. A
        // transcript with a family but no copy tag is still clustered (by family) and can carry PSVs,
        // so always show family_id when present; copy_id only when present.
        match (tx.vg_family_id, tx.vg_copy_id) {
            (Some(f), Some(c)) => {
                writeln!(fa, ">{} family_id=FAM_{} copy_id={} n_psv={}", tid, f, c, n_psv)?
            }
            (Some(f), None) => writeln!(fa, ">{} family_id=FAM_{} n_psv={}", tid, f, n_psv)?,
            _ => writeln!(fa, ">{} n_psv={}", tid, n_psv)?,
        }
        write_wrapped(&mut fa, &seq, 60)?;
        n_records += 1;

        // The emitted FASTA is in transcription orientation (reverse-complemented for `-` strand),
        // and transcript_pos indexes THAT sequence. PSV bases come from the genome-forward diff, so
        // for a `-` transcript we complement them to match the base actually at transcript_pos in the
        // emitted sequence (the user can verify base == FASTA[transcript_pos]).
        let minus = tx.strand == '-';
        for (site, tpos) in matched {
            let this_b = if minus { complement(site.this_base) } else { site.this_base };
            let siblings: String = if site.sibling_bases.is_empty() {
                ".".to_string()
            } else {
                site.sibling_bases
                    .iter()
                    .map(|&b| ((if minus { complement(b) } else { b }) as char).to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            writeln!(
                tsv,
                "{}\t{}\t{}\t{}\t{}\t{}\tFAM_{}\t{}",
                tid,
                tx.chrom,
                site.ref_pos + 1, // 1-based genomic
                tpos + 1,         // 1-based transcript
                this_b as char,
                siblings,
                site.family_id,
                site.copy_id,
            )?;
            n_rows += 1;
        }
    }
    fa.flush()?;
    tsv.flush()?;
    eprintln!(
        "[VG] PSV-FASTA: {} transcript record(s) -> {} ; {} PSV row(s) -> {}",
        n_records,
        fa_path.display(),
        n_rows,
        tsv_path.display()
    );
    Ok((n_records, n_rows))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a single-contig GenomeIndex "chrT" with the given forward sequence.
    fn make_genome(seq: &[u8]) -> GenomeIndex {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("t.fa");
        let mut content = b">chrT\n".to_vec();
        content.extend_from_slice(seq);
        content.push(b'\n');
        std::fs::write(&fa, &content).unwrap();
        GenomeIndex::from_fasta(fa.to_str().unwrap()).unwrap()
    }

    #[test]
    fn spliced_sequence_concatenates_exons_plus_strand() {
        // chrT = "AAACCCGGGTTT" ; exons [0,3) and [6,9) -> "AAA" + "GGG".
        let g = make_genome(b"AAACCCGGGTTT");
        let seq = spliced_sequence(&g, "chrT", '+', &[(0, 3), (6, 9)]);
        assert_eq!(seq, b"AAAGGG");
    }

    #[test]
    fn spliced_sequence_minus_strand_is_reverse_complement() {
        let g = make_genome(b"AAACCCGGGTTT");
        let plus = spliced_sequence(&g, "chrT", '+', &[(0, 3), (6, 9)]); // AAAGGG
        let minus = spliced_sequence(&g, "chrT", '-', &[(0, 3), (6, 9)]);
        assert_eq!(minus, crate::vg::reverse_complement(&plus));
        assert_eq!(minus, b"CCCTTT"); // revcomp(AAAGGG)
    }

    #[test]
    fn genomic_to_transcript_pos_maps_exons_and_drops_introns() {
        // exons [0,3) and [6,9): transcript = positions 0,1,2 (g0,1,2), 3,4,5 (g6,7,8).
        let exons = [(0u64, 3u64), (6u64, 9u64)];
        assert_eq!(genomic_to_transcript_pos(&exons, '+', 0), Some(0));
        assert_eq!(genomic_to_transcript_pos(&exons, '+', 2), Some(2));
        assert_eq!(genomic_to_transcript_pos(&exons, '+', 6), Some(3)); // first base of exon2
        assert_eq!(genomic_to_transcript_pos(&exons, '+', 8), Some(5));
        assert_eq!(genomic_to_transcript_pos(&exons, '+', 4), None, "intronic -> None");
        assert_eq!(genomic_to_transcript_pos(&exons, '+', 9), None, "past end -> None");
        // minus strand: emitted seq is RC, so transcript_pos = total-1 - fwd. total=6.
        assert_eq!(genomic_to_transcript_pos(&exons, '-', 0), Some(5)); // 6-1-0
        assert_eq!(genomic_to_transcript_pos(&exons, '-', 8), Some(0)); // 6-1-5
        assert_eq!(genomic_to_transcript_pos(&exons, '-', 4), None);
    }

    #[test]
    fn complement_round_trips() {
        for &b in b"ACGTacgt" {
            assert_eq!(complement(complement(b)), b);
        }
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn build_psv_map_finds_sites_from_emitted_transcripts() {
        use crate::path_extract::Transcript;
        // Two same-strand paralog copies on chrT: copyA exon [100,160), copyB exon [1000,1060),
        // identical 60bp except SNPs at exon-offset 10 and 40.
        let mut seq_a: Vec<u8> =
            b"GTCAGTTACGGATCCAAGTGCTAATCGCTGACATGGCCATTGAACCGTTAGCATGCAATG".to_vec();
        assert_eq!(seq_a.len(), 60);
        seq_a[10] = b'C';
        seq_a[40] = b'G';
        let mut seq_b = seq_a.clone();
        seq_b[10] = b'A'; // PSV 1
        seq_b[40] = b'T'; // PSV 2
        let mut g = vec![b'A'; 1200];
        g[100..160].copy_from_slice(&seq_a);
        g[1000..1060].copy_from_slice(&seq_b);
        let genome = make_genome(&g);

        let txa = Transcript {
            chrom: "chrT".into(), strand: '+', exons: vec![(100, 160)],
            vg_family_id: Some(7), vg_copy_id: Some(0), ..Default::default()
        };
        let txb = Transcript {
            chrom: "chrT".into(), strand: '+', exons: vec![(1000, 1060)],
            vg_family_id: Some(7), vg_copy_id: Some(1), ..Default::default()
        };

        let map = build_psv_map(&[txa, txb], &genome);
        let sites = map.get("chrT").expect("chrT has PSV sites");
        assert_eq!(sites.len(), 4, "2 SNPs across 2 copies → 4 site records");
        // copyA SNP at genomic 110: this_base C (seq_a[10]), sibling A (seq_b[10]).
        let a10 = sites.iter().find(|s| s.ref_pos == 110).expect("copyA site at 110");
        assert_eq!(a10.this_base, b'C');
        assert_eq!(a10.sibling_bases, vec![b'A']);
        // copyB SNP at genomic 1040: this_base T (seq_b[40]), sibling G (seq_a[40]).
        let b40 = sites.iter().find(|s| s.ref_pos == 1040).expect("copyB site at 1040");
        assert_eq!(b40.this_base, b'T');
        assert_eq!(b40.sibling_bases, vec![b'G']);
    }

    #[test]
    fn build_psv_map_skips_single_copy_family() {
        // One copy → no sibling → no PSVs.
        use crate::path_extract::Transcript;
        let g = vec![b'A'; 400];
        let genome = make_genome(&g);
        let tx = Transcript {
            chrom: "chrT".into(), strand: '+', exons: vec![(100, 160)],
            vg_family_id: Some(1), vg_copy_id: Some(0), ..Default::default()
        };
        assert!(build_psv_map(&[tx], &genome).is_empty());
    }
}
