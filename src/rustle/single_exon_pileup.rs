//! Standalone polyA-driven single-exon transcript detector.
//!
//! Operates on a bundle's read list without depending on the splice graph or
//! flow decomposition. Designed for cases where 23% of chr19 GGO missed refs
//! are single-exon transcripts (median cov=10) that the graph/flow pipeline
//! structurally cannot emit.
//!
//! Primary discrimination signal: **unspliced-read dominance** at the cluster
//! locus — for a real single-exon transcript, the bulk of supporting reads
//! cluster as unspliced and there is no competing multi-exon signal of
//! comparable depth. For noise loci (pre-mRNA fragments, retained-intron
//! artifacts, partial reads from a spliced isoform), unspliced reads are
//! a minority.
//!
//! Gated entirely by `RUSTLE_SE_PILEUP=1`. Tunable via:
//!   - `RUSTLE_SE_PILEUP_MIN_READS` (default 5) — min cluster size
//!   - `RUSTLE_SE_PILEUP_MIN_UNSPL_FRAC` (default 0.5) — cluster reads must be
//!      at least this fraction of all reads spanning the cluster region
//!   - `RUSTLE_SE_PILEUP_MIN_LEN` (default 200)
//!   - `RUSTLE_SE_PILEUP_ENDPOINT_TOL` (default 1500)
//!
//! PolyA tail detection (`has_last_exon_polya` / `has_first_exon_polyt`) was
//! evaluated as a primary discriminator but on the GGO IsoSeq BAM ~92% of
//! bundles have zero polyA-attested reads (likely pre-trimmed), so it cannot
//! drive emission. It is used only for strand assignment when present.

use crate::path_extract::Transcript;
use crate::types::{Bundle, BundleRead};

struct Config {
    min_reads: usize,
    min_unspl_frac: f64,
    min_length: u64,
    endpoint_tol: u64,
}

impl Config {
    fn from_env() -> Self {
        Self {
            min_reads: env_u("RUSTLE_SE_PILEUP_MIN_READS", 5),
            min_unspl_frac: env_f("RUSTLE_SE_PILEUP_MIN_UNSPL_FRAC", 0.5),
            min_length: env_u64("RUSTLE_SE_PILEUP_MIN_LEN", 200),
            endpoint_tol: env_u64("RUSTLE_SE_PILEUP_ENDPOINT_TOL", 1500),
        }
    }
}

fn env_u(k: &str, d: usize) -> usize {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d)
}
fn env_f(k: &str, d: f64) -> f64 {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d)
}
fn env_u64(k: &str, d: u64) -> u64 {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d)
}

#[inline]
fn abs_diff(a: u64, b: u64) -> u64 {
    if a > b { a - b } else { b - a }
}

/// Run the pileup-based SE detector on one bundle.
///
/// Returns single-exon `Transcript`s synthesized from polyA-attested
/// unspliced read clusters that don't overlap any existing multi-exon
/// transcript on the same strand.
pub fn detect(
    bundle: &Bundle,
    existing: &[Transcript],
    long_reads: bool,
) -> Vec<Transcript> {
    if std::env::var_os("RUSTLE_SE_PILEUP").is_none() {
        return Vec::new();
    }
    let cfg = Config::from_env();
    let trace = std::env::var_os("RUSTLE_SE_PILEUP_TRACE").is_some();
    if trace {
        let mono_total: usize = bundle.reads.iter()
            .filter(|r| r.junctions.is_empty() && r.exons.len() == 1).count();
        let polya_total: usize = bundle.reads.iter()
            .filter(|r| r.has_last_exon_polya || r.has_poly_end_aligned ||
                        r.has_poly_end_unaligned || r.has_first_exon_polyt ||
                        r.has_poly_start_aligned || r.has_poly_start_unaligned)
            .count();
        eprintln!(
            "[SE_PILEUP] enter bundle={}:{}-{}({}) reads={} mono={} any_polya={}",
            bundle.chrom,
            bundle.start,
            bundle.end,
            bundle.strand,
            bundle.reads.len(),
            mono_total,
            polya_total,
        );
    }

    let mono: Vec<&BundleRead> = bundle
        .reads
        .iter()
        .filter(|r| r.junctions.is_empty() && r.exons.len() == 1)
        .collect();
    if mono.len() < cfg.min_reads {
        return Vec::new();
    }

    let mut sorted = mono.clone();
    sorted.sort_by_key(|r| (r.ref_start, r.ref_end));

    let mut clusters: Vec<Vec<&BundleRead>> = Vec::new();
    let mut cur: Vec<&BundleRead> = Vec::new();
    let mut cur_anchor_start: u64 = 0;
    let mut cur_anchor_end: u64 = 0;
    for r in sorted {
        if cur.is_empty() {
            cur_anchor_start = r.ref_start;
            cur_anchor_end = r.ref_end;
            cur.push(r);
        } else {
            let ds = abs_diff(r.ref_start, cur_anchor_start);
            let de = abs_diff(r.ref_end, cur_anchor_end);
            if ds <= cfg.endpoint_tol && de <= cfg.endpoint_tol {
                cur.push(r);
            } else {
                clusters.push(std::mem::take(&mut cur));
                cur_anchor_start = r.ref_start;
                cur_anchor_end = r.ref_end;
                cur.push(r);
            }
        }
    }
    if !cur.is_empty() {
        clusters.push(cur);
    }

    let mut out = Vec::new();
    for cluster in clusters {
        let n = cluster.len();
        if n < cfg.min_reads {
            continue;
        }
        let start = cluster.iter().map(|r| r.ref_start).min().unwrap();
        let end = cluster.iter().map(|r| r.ref_end).max().unwrap();
        if end <= start || end - start < cfg.min_length {
            continue;
        }

        // Discrimination: cluster unspliced reads must be the dominant signal
        // at the locus. Count all bundle reads that overlap [start,end];
        // require cluster.len()/overlapping >= min_unspl_frac.
        let total_overlap: usize = bundle.reads.iter()
            .filter(|r| r.ref_start < end && start < r.ref_end)
            .count();
        let unspl_frac = if total_overlap > 0 { n as f64 / total_overlap as f64 } else { 0.0 };
        if unspl_frac < cfg.min_unspl_frac {
            if trace {
                eprintln!(
                    "[SE_PILEUP] skip(unspl_frac<{:.2}) {}:{}-{} n={} total_overlap={} frac={:.3}",
                    cfg.min_unspl_frac, bundle.chrom, start, end, n, total_overlap, unspl_frac
                );
            }
            continue;
        }

        // Strand: prefer polyA evidence when available, otherwise use the
        // most common read strand assignment from XS tags.
        let polya_plus = cluster.iter().filter(|r|
            r.has_last_exon_polya || r.has_poly_end_aligned || r.has_poly_end_unaligned
        ).count();
        let polya_minus = cluster.iter().filter(|r|
            r.has_first_exon_polyt || r.has_poly_start_aligned || r.has_poly_start_unaligned
        ).count();
        let strand = if polya_plus > polya_minus {
            '+'
        } else if polya_minus > polya_plus {
            '-'
        } else {
            // No polyA evidence — fall back to bundle strand if defined,
            // else the most common read XS-strand.
            let by_bundle = bundle.strand;
            if by_bundle == '+' || by_bundle == '-' {
                by_bundle
            } else {
                let p = cluster.iter().filter(|r| r.strand == '+').count();
                let m = cluster.iter().filter(|r| r.strand == '-').count();
                if p > m { '+' }
                else if m > p { '-' }
                else {
                    if trace {
                        eprintln!(
                            "[SE_PILEUP] skip(no_strand_signal) {}:{}-{} n={}",
                            bundle.chrom, start, end, n
                        );
                    }
                    continue;
                }
            }
        };

        // Reject same-strand exonic overlap with existing transcripts.
        let overlaps = existing.iter().any(|tx| {
            tx.strand == strand
                && tx.chrom == bundle.chrom
                && tx.exons.iter().any(|&(s, e)| start < e && s < end)
        });
        if overlaps {
            if trace {
                eprintln!(
                    "[SE_PILEUP] skip(overlap) {}:{}-{}({}) n={}",
                    bundle.chrom, start, end, strand, n
                );
            }
            continue;
        }

        if trace {
            eprintln!(
                "[SE_PILEUP] EMIT {}:{}-{}({}) n={} unspl_frac={:.2} polyA={}/{}",
                bundle.chrom, start, end, strand, n, unspl_frac, polya_plus, polya_minus
            );
        }

        let cov = n as f64;
        out.push(Transcript {
            chrom: bundle.chrom.clone(),
            strand,
            exons: vec![(start, end)],
            coverage: cov,
            exon_cov: vec![cov],
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("se_pileup".to_string()),
            is_longread: long_reads,
            longcov: cov,
            bpcov_cov: 0.0,
            all_strand_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
            alt_tts_end: false,
            vg_family_id: None,
            vg_copy_id: None,
            vg_family_size: None,
            intron_low: Vec::new(),
            synthetic: false,
            rescue_class: None,
            raw_flow_sum: 0.0,

});
    }
    out
}
