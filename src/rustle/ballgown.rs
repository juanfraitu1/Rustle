//! Ballgown-style tabular output (t_data.ctab, e_data.ctab) for downstream analysis.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::path_extract::Transcript;

/// Write Ballgown-style tables into `dir`: t_data.ctab (transcript-level), e_data.ctab (exon-level).
/// Coordinates are 1-based. Creates directory if needed.
pub fn write_ballgown<P: AsRef<Path>>(
    transcripts: &[Transcript],
    dir: P,
    label: &str,
) -> std::io::Result<()> {
    let dir = dir.as_ref();
    std::fs::create_dir_all(dir)?;

    // t_data.ctab: t_id, gene_id, chr, strand, start, end, num_exons, length, coverage, tpm, fpkm
    let t_path = dir.join("t_data.ctab");
    let mut tw = BufWriter::new(File::create(t_path)?);
    writeln!(
        tw,
        "t_id\tgene_id\tchr\tstrand\tstart\tend\tnum_exons\tlength\tcoverage\ttpm\tfpkm"
    )?;
    for (i, tx) in transcripts.iter().enumerate() {
        let t_id = format!("{}.{}", label, i + 1);
        let start = tx.exons.first().map(|(s, _)| s + 1).unwrap_or(0);
        let end = tx.exons.last().map(|(_, e)| *e).unwrap_or(0);
        let num_exons = tx.exons.len();
        let length: u64 = tx.exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
        writeln!(
            tw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
            t_id,
            t_id,
            tx.chrom,
            tx.strand,
            start,
            end,
            num_exons,
            length,
            tx.coverage,
            tx.tpm,
            tx.fpkm
        )?;
    }
    tw.flush()?;

    // e_data.ctab: e_id, chr, strand, start, end, t_id, gene_id, exon_number, length, coverage
    let e_path = dir.join("e_data.ctab");
    let mut ew = BufWriter::new(File::create(e_path)?);
    writeln!(
        ew,
        "e_id\tchr\tstrand\tstart\tend\tt_id\tgene_id\texon_number\tlength\tcoverage"
    )?;
    for (i, tx) in transcripts.iter().enumerate() {
        let t_id = format!("{}.{}", label, i + 1);
        for (j, (s, e)) in tx.exons.iter().enumerate() {
            let e_id = format!("{}.{}", t_id, j + 1);
            let length = e.saturating_sub(*s);
            writeln!(
                ew,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                e_id,
                tx.chrom,
                tx.strand,
                s + 1,
                e,
                t_id,
                t_id,
                j + 1,
                length,
                tx.coverage
            )?;
        }
    }
    ew.flush()?;

    // i_data.ctab: introns between consecutive exons (1-based coordinates)
    // i_id, chr, strand, start, end, t_id, gene_id, exon1_number, exon2_number, length, coverage
    let i_path = dir.join("i_data.ctab");
    let mut iw = BufWriter::new(File::create(i_path)?);
    writeln!(
        iw,
        "i_id\tchr\tstrand\tstart\tend\tt_id\tgene_id\texon1_number\texon2_number\tlength\tcoverage"
    )?;
    for (i, tx) in transcripts.iter().enumerate() {
        let t_id = format!("{}.{}", label, i + 1);
        for j in 0..tx.exons.len().saturating_sub(1) {
            let (_, e1_end) = tx.exons[j];
            let (e2_start, _) = tx.exons[j + 1];
            let intron_start = e1_end; // 0-based first intron base
            let intron_end = e2_start; // 0-based last intron base (exclusive in typical 0-based)
            let length = intron_end.saturating_sub(intron_start);
            let i_id = format!("{}.i{}", t_id, j + 1);
            writeln!(
                iw,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                i_id,
                tx.chrom,
                tx.strand,
                intron_start + 1,
                intron_end,
                t_id,
                t_id,
                j + 1,
                j + 2,
                length,
                tx.coverage
            )?;
        }
    }
    iw.flush()?;

    Ok(())
}
