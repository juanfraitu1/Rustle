//! Output splice junctions from final transcripts (for visualization / downstream).

use crate::types::DetHashSet as HashSet;
use std::io::Write;

use crate::path_extract::Transcript;

/// One splice junction: donor (exon end) and acceptor (exon start). 0-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct OutputJunction {
    pub donor: u64,
    pub acceptor: u64,
}

/// Collect unique (chr, strand, donor, acceptor) from all transcripts.
pub fn junctions_from_transcripts(transcripts: &[Transcript]) -> Vec<(String, char, OutputJunction)> {
    let mut seen: HashSet<(String, char, u64, u64)> = Default::default();
    let mut out: Vec<(String, char, OutputJunction)> = Vec::new();
    for tx in transcripts {
        for j in 0..tx.exons.len().saturating_sub(1) {
            let (_, e1_end) = tx.exons[j];
            let (e2_start, _) = tx.exons[j + 1];
            let key = (tx.chrom.clone(), tx.strand, e1_end, e2_start);
            if seen.insert((tx.chrom.clone(), tx.strand, e1_end, e2_start)) {
                out.push((
                    tx.chrom.clone(),
                    tx.strand,
                    OutputJunction {
                        donor: e1_end,
                        acceptor: e2_start,
                    },
                ));
            }
        }
    }
    out
}

/// Write junctions to a tab-separated file: chr, strand, donor, acceptor (1-based for display).
pub fn write_junctions_table<W: Write>(
    transcripts: &[Transcript],
    writer: &mut W,
) -> std::io::Result<()> {
    writeln!(writer, "chr\tstrand\tdonor\tacceptor")?;
    let juncs = junctions_from_transcripts(transcripts);
    for (chr, strand, j) in juncs {
        writeln!(writer, "{}\t{}\t{}\t{}", chr, strand, j.donor + 1, j.acceptor + 1)?;
    }
    Ok(())
}
