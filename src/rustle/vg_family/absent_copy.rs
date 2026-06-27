//! Admission gate for reference-ABSENT (collapsed) copies.
//!
//! A `CollapsedCandidate` (from `copy_split`) passes through five ordered gates and either
//! becomes an admitted `Copy` (a synthetic `DenovoTranscript` carrying the copy's distinguishing
//! alleles) or a `DnaNeeds` record flagged for downstream DNA-level validation.  All five gates
//! must pass — first failure determines the reason string.
//!
//! Gate order:
//!   1. `n_clusters >= min_clusters` — enough co-varying clusters for a multi-copy call.
//!   2. `min_p_distinct` — the candidate's allele vector is certifiably distinct from the host.
//!   3. `strand_symmetric_spectrum` — NOT a pure A→G editing artefact.
//!   4. `collapsed_copy_to_transcript_from_host_seq` — alleles are placeable (substitution-only).
//!   5. Genome remap identity `< remap_max_identity` — the synthetic copy does not remap
//!      back to its host locus with ≥ 98 % identity (i.e. it is genuinely divergent).
//!
//! The hermetic `admit_candidate_with_remap` takes the already-spliced host transcript and an
//! injected `remap_identity` closure so tests never shell minimap2.  The real entry point
//! `admit_candidate` rebuilds the host spliced sequence from the genome (using the copy's intron
//! chain), then delegates to the hermetic variant with a real minimap2 closure.

use std::collections::BTreeMap;

use crate::genome::GenomeIndex;
use crate::vg_family::copy_assign_pipeline::gen2off;
use crate::vg_family::copy_split::{
    collapsed_copy_to_transcript_from_host_seq, min_p_distinct, strand_symmetric_spectrum,
    CollapsedCandidate,
};
use crate::vg_family::family_detect::DenovoTranscript;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Tunable thresholds for the admission gate.
pub struct AbsentCopyParams {
    /// Per-base sequencing error rate used in `min_p_distinct` (default 0.003 = HiFi).
    pub error_rate: f64,
    /// Bonferroni-corrected significance threshold for `min_p_distinct` (default 1e-3).
    pub alpha: f64,
    /// Minimum number of identifiable co-varying clusters at the locus (default 3).
    pub min_clusters: usize,
    /// Upper bound on remap identity; above this threshold the copy maps back too cleanly
    /// to be a genuine divergent paralog (default 0.98 = 98 %).
    pub remap_max_identity: f64,
}

impl Default for AbsentCopyParams {
    fn default() -> Self {
        Self {
            error_rate: 0.003,
            alpha: 1e-3,
            min_clusters: 3,
            remap_max_identity: 0.98,
        }
    }
}

/// A candidate that failed one or more gates and requires DNA-level validation.
#[derive(Clone, Debug)]
pub struct DnaNeedsRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_clusters: usize,
    /// Human-readable string identifying the first gate that fired.
    pub reason: String,
    pub read_count: usize,
}

/// Result of the admission gate.
pub enum Admission {
    /// All gates passed — the `DenovoTranscript` is the synthetic collapsed copy.
    Copy(DenovoTranscript),
    /// At least one gate failed — the record is flagged for DNA validation.
    DnaNeeds(DnaNeedsRecord),
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn dna_needs(c: &CollapsedCandidate, reason: &str) -> Admission {
    Admission::DnaNeeds(DnaNeedsRecord {
        chrom: c.chrom.clone(),
        start: c.start,
        end: c.end,
        n_clusters: c.n_clusters,
        reason: reason.into(),
        read_count: c.iso.read_count,
    })
}

/// Derive the host's own bases at every PSV position from the already-spliced `host_spliced`
/// transcript (whose `gen2off` map covers exactly the exons built with the copy's intron chain).
fn host_alleles_at(host_spliced: &DenovoTranscript, psv_pos: &[u64]) -> Vec<Option<u8>> {
    let g2o: BTreeMap<u64, usize> = gen2off(host_spliced);
    psv_pos
        .iter()
        .map(|&pos| g2o.get(&pos).and_then(|&off| host_spliced.seq.get(off).copied()))
        .collect()
}

// ---------------------------------------------------------------------------
// Core gate logic (hermetic — remap is injected)
// ---------------------------------------------------------------------------

/// Hermetic admission gate.  `host_spliced` must already carry the COPY's intron chain and
/// the corresponding spliced sequence (so `gen2off` agrees with `host_spliced.seq`).
/// `remap_identity` is a closure `|seq: &[u8]| -> Option<f64>` that returns the best
/// alignment identity against the genome (or `None` if there is no alignment).  In tests,
/// pass `|_seq| Some(0.5)` to bypass minimap2 entirely.
pub(crate) fn admit_candidate_with_remap<F>(
    cand: &CollapsedCandidate,
    host_spliced: &DenovoTranscript,
    p: &AbsentCopyParams,
    remap_identity: F,
) -> Admission
where
    F: Fn(&[u8]) -> Option<f64>,
{
    // Gate 1: cluster count.
    if cand.n_clusters < p.min_clusters {
        return dna_needs(
            cand,
            &format!("<{} clusters (copy-vs-allele needs DNA)", p.min_clusters),
        );
    }

    // Gate 2: min_p distinct from host.
    let host_alleles = host_alleles_at(host_spliced, &cand.psv_pos);
    if !min_p_distinct(&cand.iso.allele_vector, &host_alleles, p.error_rate, p.alpha) {
        return dna_needs(cand, "not min_p-distinct from host");
    }

    // Gate 3: strand-symmetric spectrum (reject pure A→G editing artefact).
    if !strand_symmetric_spectrum(&host_alleles, &cand.iso.allele_vector) {
        return dna_needs(cand, "pure A->G spectrum (editing-like)");
    }

    // Gate 4: build the overlay transcript.
    let t = match collapsed_copy_to_transcript_from_host_seq(
        &cand.iso,
        &cand.psv_pos,
        host_spliced,
    ) {
        Some(t) => t,
        None => return dna_needs(cand, "consensus unplaceable (indel/intron)"),
    };

    // Gate 5: genome remap identity.
    match remap_identity(&t.seq) {
        Some(id) if id < p.remap_max_identity => Admission::Copy(t),
        Some(_) => dna_needs(cand, ">=98% remap identity (paralog-leak or het)"),
        None => dna_needs(cand, "no homology on remap"),
    }
}

// ---------------------------------------------------------------------------
// Real entry point
// ---------------------------------------------------------------------------

/// Admission gate with real genome-based host-seq reconstruction and minimap2 remap.
///
/// `fasta_path` is the genome FASTA file (required by minimap2; `GenomeIndex` does not expose
/// its path so it is threaded in separately).  If the host spliced sequence cannot be built
/// (e.g. non-canonical junction in the copy's intron chain), routes to `DnaNeeds`.
pub fn admit_candidate(
    cand: &CollapsedCandidate,
    host: &DenovoTranscript,
    genome: &GenomeIndex,
    fasta_path: &str,
    p: &AbsentCopyParams,
) -> Admission {
    use crate::vg_family::denovo_assemble::build_spliced_seq;

    // Rebuild host spliced seq using the COPY's intron chain (not the host's own chain), so that
    // `gen2off` inside the hermetic gate agrees with the seq bytes.
    let (seq, strand) = match build_spliced_seq(
        genome,
        &host.chrom,
        host.start,
        host.end,
        &cand.iso.intron_chain,
    ) {
        Some(v) => v,
        None => return dna_needs(cand, "host sequence unbuildable"),
    };
    let host_spliced = DenovoTranscript {
        seq,
        strand,
        introns: cand.iso.intron_chain.clone(),
        ..host.clone()
    };

    let fp = fasta_path.to_owned();
    admit_candidate_with_remap(cand, &host_spliced, p, move |seq| {
        remap_identity_minimap2(seq, &fp)
    })
}

// ---------------------------------------------------------------------------
// minimap2 remap helper
// ---------------------------------------------------------------------------

/// Shell minimap2 (`-cx asm20 --eqx`) to remap `query_seq` against `fasta_path`.
/// Returns `Some(identity)` for the best alignment (by alignment block length), `None` if there
/// is no alignment.  Identity is derived from the `de:f:<x>` divergence tag when present
/// (`identity = 1 − x`), falling back to `matches / aln_block_len` (PAF cols 10/11).
fn remap_identity_minimap2(query_seq: &[u8], fasta_path: &str) -> Option<f64> {
    use std::io::Write;

    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    // Nonce derived from content length so concurrent calls with different queries don't collide.
    let nonce = query_seq.len().wrapping_mul(1_000_003usize);
    let qpath = dir.join(format!("rustle_absent_q_{pid}_{nonce}.fa"));

    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let _cl = Cleanup(qpath.clone());

    {
        let mut q = std::fs::File::create(&qpath).ok()?;
        q.write_all(b">q\n").ok()?;
        q.write_all(query_seq).ok()?;
        q.write_all(b"\n").ok()?;
    }

    let out = std::process::Command::new(&mm2)
        .args(["-cx", "asm20", "--eqx", "--secondary=no", "-t", "1"])
        .arg(fasta_path)
        .arg(&qpath)
        .output()
        .ok()?;

    if !out.status.success() {
        return None;
    }

    let text = String::from_utf8_lossy(&out.stdout);
    let mut best_id: Option<f64> = None;
    let mut best_span: u64 = 0;

    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        // PAF col 10 = number of matching bases; col 11 = alignment block length.
        let aln_block_len: u64 = match f[10].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if aln_block_len == 0 {
            continue;
        }
        // Prefer the de:f: divergence tag (more accurate for asm20 output).
        let id = if let Some(de_str) = f.iter().find_map(|x| x.strip_prefix("de:f:")) {
            let div: f64 = de_str.parse().unwrap_or(1.0);
            1.0 - div
        } else {
            let matches: u64 = f[9].parse().unwrap_or(0);
            matches as f64 / aln_block_len as f64
        };

        if aln_block_len > best_span {
            best_span = aln_block_len;
            best_id = Some(id);
        }
    }

    best_id
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::copy_split::{CollapsedCandidate, CopyIsoform};
    use crate::vg_family::family_detect::DenovoTranscript;

    fn host() -> DenovoTranscript {
        DenovoTranscript {
            tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            n_reads: 9,
            strand: '+',
            introns: vec![],
            seq: b"AAAAAAAAAA".to_vec(),
        }
    }

    /// Gate 1 fires: only 2 clusters (< min_clusters=3).  The remap closure returns 0.5
    /// (would otherwise be divergent enough to pass gate 5), so the only failure is gate 1.
    #[test]
    fn admit_rejects_two_cluster_substitution_only_as_dna_needs() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'C'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 2, // only 2 clusters
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("clusters"),
                "expected 'clusters' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// All five gates pass: 3 clusters, allele vector is min_p-distinct from host (3 differing
    /// columns → p = (ε/3)³ ≪ α), mixed spectrum (C/T/C not A→G), overlay is buildable
    /// (single-exon identity map), remap returns 0.5 < 0.98.
    #[test]
    fn admit_accepts_three_cluster_distinct_divergent() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::Copy(t) => assert_eq!(t.chrom, "c1"),
            _ => panic!("expected Copy"),
        }
    }

    /// Gate 2 fires: allele vector identical to host at all observed columns → not min_p-distinct.
    #[test]
    fn admit_rejects_non_distinct_from_host() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'A'), Some(b'A'), Some(b'A')], // same as host
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("min_p-distinct"),
                "expected 'min_p-distinct' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// Gate 3 fires: pure A→G spectrum (RNA-editing artefact).
    #[test]
    fn admit_rejects_pure_ag_spectrum() {
        // Host is all A; candidate is all G at those positions → pure A→G.
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'G'), Some(b'G'), Some(b'G')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("A->G"),
                "expected 'A->G' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// Gate 5 fires: remap identity ≥ remap_max_identity (0.99 ≥ 0.98).
    #[test]
    fn admit_rejects_high_remap_identity() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.99));
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("remap identity"),
                "expected 'remap identity' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// Gate 5 fires: no alignment returned (None).
    #[test]
    fn admit_rejects_no_remap_alignment() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "c1".into(),
            start: 100,
            end: 110,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5,
                identifiable: true,
            },
            psv_pos: vec![102, 104, 106],
            n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| None);
        match got {
            Admission::DnaNeeds(r) => assert!(
                r.reason.contains("no homology"),
                "expected 'no homology' in reason, got: {}",
                r.reason
            ),
            _ => panic!("expected DnaNeeds"),
        }
    }

    /// DnaNeedsRecord carries correct metadata from the candidate.
    #[test]
    fn dna_needs_record_fields_from_candidate() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(),
            chrom: "chr7".into(),
            start: 500,
            end: 600,
            iso: CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C')],
                read_count: 12,
                identifiable: true,
            },
            psv_pos: vec![502],
            n_clusters: 1, // gate 1 fires
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_| Some(0.5));
        match got {
            Admission::DnaNeeds(r) => {
                assert_eq!(r.chrom, "chr7");
                assert_eq!(r.start, 500);
                assert_eq!(r.end, 600);
                assert_eq!(r.n_clusters, 1);
                assert_eq!(r.read_count, 12);
            }
            _ => panic!("expected DnaNeeds"),
        }
    }
}
