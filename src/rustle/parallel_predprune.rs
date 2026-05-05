//! Parallel prediction-pruning pass — default off.
//!
//! Targets two parity gaps observed on GGO_19 vs StringTie:
//!   - **m-class** (24 / ~1900 tx): rustle's flow extracts BOTH the spliced
//!     and a merged-exon path; ST only emits the spliced path. The merged
//!     path is a retained-intron variant that ST's flow doesn't enumerate.
//!   - **k-class** (26 tx): rustle extends transcripts past ST's reference
//!     through low-support terminal junctions, typically driven by a single
//!     long-spanning read.
//!
//! ## Hypothesis
//!
//! ST applies "junction-support relative to gene's max coverage" gating
//! at terminal regions, and rejects RI variants where the body coverage
//! is small relative to the corresponding spliced junction's support.
//!
//! Both signals can be computed from the post-extraction transcript set
//! alone, without re-running the flow algorithm:
//!   - per-locus `gene_max_cov` = max transcript cov in the cluster
//!   - per-junction `support` = sum of `cov × len` across transcripts that
//!     contain the junction (proxy for its read backing)
//!   - per-merged-exon body coverage = the cov of the variant containing
//!     the merged exon (assumes flow has correctly redistributed reads)
//!
//! ## Approach
//!
//! Subtractive — trim or drop transcripts; never add. Two passes:
//!   1. `m_kill_pass`: for each tx pair (n1, n2) at the same locus where
//!      n2's intron-set ⊂ n1's AND n2 has an exon strictly containing
//!      n1's missing introns (= classic m-class signature), kill n2 only
//!      if its cov is below `m_kill_frac × n1.cov`. This is the same
//!      cov-ratio rule explored in path_extract, but here it runs over
//!      finalized transcripts so we can use the actual emitted cov rather
//!      than mid-flow estimates.
//!   2. `k_trim_pass`: for each tx, walk terminal junctions inward from
//!      the ends. If both flanking exons have per-base cov < `k_min_support`
//!      (absolute threshold, default 3.0 reads/base — a "supported by ≥3
//!      reads" gate), truncate the tx at that junction. Absolute threshold
//!      avoids the cluster-dominant-isoform pitfall: a high-cov sibling
//!      doesn't push the threshold up and accidentally clip neighbors with
//!      legitimate moderate-cov UTRs.
//!
//! ## Gating
//!
//! - `RUSTLE_PARALLEL_PRUNE=1` enables the pass entirely.
//! - `RUSTLE_PARALLEL_PRUNE_M_FRAC` cov-ratio for m-kill (default 0.30).
//! - `RUSTLE_PARALLEL_PRUNE_K_MIN_SUPPORT` per-base cov threshold for
//!   k-trim (default 3.0).
//! - `RUSTLE_PARALLEL_PRUNE_K_MAX_LEN` only trim terminal exons whose
//!   length is ≤ this (default 200 bp); avoids killing legitimate UTRs.
//! - `RUSTLE_PARALLEL_PRUNE_DRY_RUN=1` logs would-fire actions without acting.
//! - `RUSTLE_PARALLEL_PRUNE_M_OFF=1` skips m-kill, runs only k-trim.
//! - `RUSTLE_PARALLEL_PRUNE_K_OFF=1` skips k-trim, runs only m-kill.

use crate::path_extract::Transcript;
use std::collections::{HashMap, HashSet};

const DEFAULT_M_FRAC: f64 = 0.30;
const DEFAULT_K_MIN_SUPPORT: f64 = 3.0;
const DEFAULT_K_MAX_LEN: u64 = 200;

/// Run both passes. Returns (filtered, n_m_killed, n_k_trimmed_first, n_k_trimmed_last).
pub fn parallel_prune(
    mut transcripts: Vec<Transcript>,
    verbose: bool,
) -> (Vec<Transcript>, usize, usize, usize) {
    let dry_run = std::env::var_os("RUSTLE_PARALLEL_PRUNE_DRY_RUN").is_some();
    let m_off = std::env::var_os("RUSTLE_PARALLEL_PRUNE_M_OFF").is_some();
    let k_off = std::env::var_os("RUSTLE_PARALLEL_PRUNE_K_OFF").is_some();

    let m_frac: f64 = std::env::var("RUSTLE_PARALLEL_PRUNE_M_FRAC")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_M_FRAC);
    let k_min_support: f64 = std::env::var("RUSTLE_PARALLEL_PRUNE_K_MIN_SUPPORT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_K_MIN_SUPPORT);
    let k_max_len: u64 = std::env::var("RUSTLE_PARALLEL_PRUNE_K_MAX_LEN")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_K_MAX_LEN);

    // Group by locus (chrom + overlapping span, same strand).
    let clusters = group_by_locus(&transcripts);

    let mut to_kill: HashSet<usize> = HashSet::new();
    let mut k_first = 0usize;
    let mut k_last = 0usize;

    for cluster in &clusters {
        if cluster.len() < 2 {
            continue;
        }
        // Pass 1: m-kill
        if !m_off {
            for &n1 in cluster {
                if to_kill.contains(&n1) {
                    continue;
                }
                let t1 = &transcripts[n1];
                let n1_introns = intron_set(t1);
                if n1_introns.is_empty() {
                    continue;
                }
                for &n2 in cluster {
                    if n1 == n2 || to_kill.contains(&n2) {
                        continue;
                    }
                    let t2 = &transcripts[n2];
                    if t2.exons.len() >= t1.exons.len() {
                        continue;
                    }
                    if t2.coverage >= m_frac * t1.coverage {
                        continue;
                    }
                    if t2.ref_transcript_id.is_some() {
                        continue;
                    }
                    let n2_introns = intron_set(t2);
                    if !n2_introns.is_subset(&n1_introns) {
                        continue;
                    }
                    let missing: Vec<(u64, u64)> =
                        n1_introns.difference(&n2_introns).copied().collect();
                    if missing.is_empty() {
                        continue;
                    }
                    let all_spanned = missing.iter().all(|&(donor, acc)| {
                        t2.exons.iter().any(|&(es, ee)| es < donor && ee > acc)
                    });
                    if !all_spanned {
                        continue;
                    }
                    if !dry_run {
                        to_kill.insert(n2);
                    }
                    if verbose {
                        eprintln!(
                            "[PARALLEL_PRUNE_M] killer={} cov={:.2} victim={} cov={:.2} missing={}",
                            tx_id(t1),
                            t1.coverage,
                            tx_id(t2),
                            t2.coverage,
                            missing.len()
                        );
                    }
                }
            }
        }

        // Pass 2: k-trim using absolute support threshold
        if !k_off {
            for &i in cluster {
                if to_kill.contains(&i) {
                    continue;
                }
                if transcripts[i].ref_transcript_id.is_some() {
                    continue;
                }
                if transcripts[i].exons.len() < 2 {
                    continue;
                }
                let (trimmed_first, trimmed_last) = k_trim_one(
                    &mut transcripts[i],
                    k_min_support,
                    k_max_len,
                    dry_run,
                    verbose,
                );
                k_first += trimmed_first;
                k_last += trimmed_last;
            }
        }
    }

    let n_killed = to_kill.len();
    if !to_kill.is_empty() {
        let kept: Vec<Transcript> = transcripts
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !to_kill.contains(i))
            .map(|(_, t)| t)
            .collect();
        transcripts = kept;
    }
    if verbose {
        eprintln!(
            "[PARALLEL_PRUNE] m_killed={} k_trim_first={} k_trim_last={} dry_run={}",
            n_killed, k_first, k_last, dry_run
        );
    }
    (transcripts, n_killed, k_first, k_last)
}

/// Iteratively peel terminal exons whose junction support (proxy: min
/// per-base cov of the two flanking exons) is below `min_support` AND
/// whose own length is ≤ `max_len`. The size gate avoids killing real
/// long UTRs that happen to be lower-cov; we target the "isolated long
/// read drove a multi-junction extension" signature where each upstream
/// exon is short (60-150 bp) and supported by 1-3 reads.
fn k_trim_one(
    tx: &mut Transcript,
    min_support: f64,
    max_len: u64,
    dry_run: bool,
    verbose: bool,
) -> (usize, usize) {
    let mut trimmed_first = 0;
    let mut trimmed_last = 0;
    let max_iters = 32usize;

    // Trim from front
    for _ in 0..max_iters {
        if tx.exons.len() < 2 {
            break;
        }
        if tx.exon_cov.len() != tx.exons.len() {
            break;
        }
        let (fs, fe) = tx.exons[0];
        let flen = fe.saturating_sub(fs);
        if flen > max_len {
            break;
        }
        // Junction support proxy: min cov of the two exons flanking
        // the first junction (= reads spanning that junction).
        let support = tx.exon_cov[0].min(tx.exon_cov[1]);
        if support >= min_support {
            break;
        }
        if verbose {
            eprintln!(
                "[PARALLEL_PRUNE_K] trim_first tx={} support={:.2} flen={} thr={:.2}",
                tx_id(tx),
                support,
                flen,
                min_support
            );
        }
        if !dry_run {
            tx.exons.remove(0);
            tx.exon_cov.remove(0);
            if !tx.intron_low.is_empty() {
                tx.intron_low.remove(0);
            }
        }
        trimmed_first += 1;
        if dry_run {
            break;
        }
    }
    // Trim from back
    for _ in 0..max_iters {
        if tx.exons.len() < 2 {
            break;
        }
        if tx.exon_cov.len() != tx.exons.len() {
            break;
        }
        let n = tx.exons.len();
        let (ls, le) = tx.exons[n - 1];
        let llen = le.saturating_sub(ls);
        if llen > max_len {
            break;
        }
        let support = tx.exon_cov[n - 1].min(tx.exon_cov[n - 2]);
        if support >= min_support {
            break;
        }
        if verbose {
            eprintln!(
                "[PARALLEL_PRUNE_K] trim_last tx={} support={:.2} llen={} thr={:.2}",
                tx_id(tx),
                support,
                llen,
                min_support
            );
        }
        if !dry_run {
            tx.exons.pop();
            tx.exon_cov.pop();
            if !tx.intron_low.is_empty() {
                tx.intron_low.pop();
            }
        }
        trimmed_last += 1;
        if dry_run {
            break;
        }
    }
    (trimmed_first, trimmed_last)
}

/// Group transcripts by chromosome + strand + overlapping span.
fn group_by_locus(transcripts: &[Transcript]) -> Vec<Vec<usize>> {
    let mut by_key: HashMap<(String, char), Vec<(usize, u64, u64)>> = HashMap::new();
    for (i, t) in transcripts.iter().enumerate() {
        let Some(s) = t.exons.first().map(|e| e.0) else {
            continue;
        };
        let Some(e) = t.exons.last().map(|e| e.1) else {
            continue;
        };
        by_key
            .entry((t.chrom.clone(), t.strand))
            .or_default()
            .push((i, s, e));
    }
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for (_, mut v) in by_key {
        v.sort_by_key(|x| (x.1, x.2));
        let mut cur: Vec<usize> = Vec::new();
        let mut cur_end = 0u64;
        for (i, s, e) in v {
            if cur.is_empty() || s >= cur_end {
                if !cur.is_empty() {
                    clusters.push(std::mem::take(&mut cur));
                }
                cur_end = e;
            } else {
                cur_end = cur_end.max(e);
            }
            cur.push(i);
        }
        if !cur.is_empty() {
            clusters.push(cur);
        }
    }
    clusters
}

fn intron_set(t: &Transcript) -> HashSet<(u64, u64)> {
    t.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

fn tx_id(t: &Transcript) -> String {
    t.transcript_id
        .clone()
        .unwrap_or_else(|| format!("{}:{}-{}", t.chrom, t.exons.first().map(|e| e.0).unwrap_or(0), t.exons.last().map(|e| e.1).unwrap_or(0)))
}
