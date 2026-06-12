//! -G2 annotation -> position-agnostic paralog families.
//! Clusters annotation transcripts into families by SHARED-EXON SEQUENCE similarity
//! (never coordinate proximity), gated by RunConfig.family_exon_similarity.

use anyhow::Result;
use std::path::Path;

/// One annotated transcript = one copy in a (potential) family.
#[derive(Clone, Debug)]
pub struct CopyStructure {
    pub copy_id: String,
    pub chrom: String,
    pub strand: char,
    /// 0-based half-open exons.
    pub exons: Vec<(u64, u64)>,
}

/// A paralog family: >=2 copies merged by shared-exon similarity, with the achieved similarity.
#[derive(Clone, Debug)]
pub struct AnnotationFamily {
    pub family_id: String,
    pub copies: Vec<CopyStructure>,
    pub achieved_min_sim: f64,
    pub achieved_mean_sim: f64,
}

/// Concatenated exon sequence of a copy, in genomic (exon-list) order.
///
/// Fetches each exon's bases from the genome using 0-based half-open coordinates
/// (matching `CopyStructure.exons`). Returns `None` if any exon is entirely out
/// of range for its chromosome. Strand is NOT considered here — all copies are
/// returned in the same forward-genomic convention so similarity is computed on
/// a consistent basis.
pub fn copy_sequence(
    copy: &CopyStructure,
    genome: &crate::genome::GenomeIndex,
) -> Option<Vec<u8>> {
    let mut seq = Vec::new();
    for &(s, e) in &copy.exons {
        let part = genome.fetch_sequence(&copy.chrom, s, e)?;
        seq.extend_from_slice(&part);
    }
    Some(seq)
}

/// Load -G2 GTF into copy structures (one per transcript).
pub fn load_copies<P: AsRef<Path>>(path: P) -> Result<Vec<CopyStructure>> {
    let refs = crate::reference_gtf::parse_reference_gtf(path)?;
    Ok(refs
        .into_iter()
        .map(|r| CopyStructure {
            copy_id: r.id,
            chrom: r.chrom,
            strand: r.strand,
            exons: r.exons,
        })
        .collect())
}
