//! Cross-strand predcluster pass — parallel module, default off.
//!
//! Runs as a post-processing step on the final `all_transcripts` set after
//! every per-strand bundle has been processed. Groups transcripts by
//! genomic span overlap (regardless of strand), then applies StringTie's
//! `retainedintron` rule with `n1` allowed to be on either strand.
//!
//! ## Why
//!
//! StringTie's bundles are cross-strand: at a given locus, both `+` and
//! `-` strand transcripts feed `pairwise_predord` together. A high-cov
//! antisense neighbor can drive `retainedintron` LEADING-return-1 kills
//! against low-cov sense-strand transcripts (e.g. RSTL.135 + strand at
//! 44669742-44715975 vs the - strand neighbor at 44646789-44670832,
//! ST cov=182.24).
//!
//! Rustle's bundles are per-strand (`graph_bundle.strand` at
//! `pipeline.rs:10736`), so the antisense `n1` candidate never enters
//! the per-strand predcluster scope. This module patches that gap.
//!
//! ## Approach
//!
//! Subtractive only — never adds transcripts. Operates on the post-pipeline
//! set, groups by overlap, and applies cross-strand `retained_intron_score`
//! with `n1` sorted by `cov` desc within each group.
//!
//! Without per-region bpcov (it's gone by this point), we approximate
//! the lowintron flag: for high-cov `n1` (cov >= MIN_LOWINTRON_COV), assume
//! all introns are "low". For lower-cov `n1`, do not run this rule (its
//! kills would be too aggressive without bpcov verification). This matches
//! ST behavior in spirit — only dominant transcripts drive cross-strand
//! kills.
//!
//! ## Gating
//!
//! - `RUSTLE_CROSS_STRAND_KRI=1` — enable the pass.
//! - `RUSTLE_CROSS_STRAND_DRY_RUN=1` — log would-fire kills without acting.
//! - `RUSTLE_CROSS_STRAND_MIN_COV=N` — override `MIN_LOWINTRON_COV` (default 50.0).

use crate::path_extract::Transcript;
use crate::transcript_filter::retained_intron_score;

const DEFAULT_MIN_LOWINTRON_COV: f64 = 50.0;

/// Group transcript indices by chromosome + overlapping span.
/// Returns Vec<Vec<usize>> where each inner Vec is one cluster.
fn group_by_overlap(transcripts: &[Transcript]) -> Vec<Vec<usize>> {
    if transcripts.is_empty() {
        return Vec::new();
    }
    // (idx, chrom, start, end)
    let mut indexed: Vec<(usize, &str, u64, u64)> = transcripts
        .iter()
        .enumerate()
        .filter_map(|(i, t)| {
            let s = t.exons.first().map(|e| e.0)?;
            let e = t.exons.last().map(|e| e.1)?;
            Some((i, t.chrom.as_str(), s, e))
        })
        .collect();
    indexed.sort_by(|a, b| {
        a.1.cmp(b.1)
            .then(a.2.cmp(&b.2))
            .then(b.3.cmp(&a.3)) // longer ends first within tie
    });

    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut cur: Vec<usize> = Vec::new();
    let mut cur_chrom: &str = "";
    let mut cur_end: u64 = 0;
    for (i, chrom, s, e) in indexed {
        if cur.is_empty() || chrom != cur_chrom || s >= cur_end {
            if !cur.is_empty() {
                clusters.push(std::mem::take(&mut cur));
            }
            cur_chrom = chrom;
            cur_end = e;
        } else {
            cur_end = cur_end.max(e);
        }
        cur.push(i);
    }
    if !cur.is_empty() {
        clusters.push(cur);
    }
    clusters
}

/// Heuristic: assume all of `n1`'s introns are "low" if `n1.cov` is at or
/// above `min_lowintron_cov`. Otherwise return all-false (no kills via
/// this n1). Avoids needing per-region bpcov which is unavailable here.
fn assumed_lowintron(n1: &Transcript, min_lowintron_cov: f64) -> Vec<bool> {
    if n1.exons.len() < 2 {
        return Vec::new();
    }
    let n_introns = n1.exons.len() - 1;
    if n1.coverage >= min_lowintron_cov {
        vec![true; n_introns]
    } else {
        vec![false; n_introns]
    }
}

/// Mirror ST's `overlaps[n1, n2]` gate: at least one exon of `a` and one
/// exon of `b` share a base. Span-only overlap is NOT enough (a transcript
/// whose last exon ends just before another's first exon shares the bundle
/// region but not actual exonic territory).
fn exons_share_a_base(a: &Transcript, b: &Transcript) -> bool {
    let mut i = 0;
    let mut j = 0;
    while i < a.exons.len() && j < b.exons.len() {
        let (as_, ae) = a.exons[i];
        let (bs, be) = b.exons[j];
        // 0-based half-open overlap: max(start) < min(end).
        if as_.max(bs) < ae.min(be) {
            return true;
        }
        // Advance whichever ends first.
        if ae <= be {
            i += 1;
        } else {
            j += 1;
        }
    }
    false
}

/// Apply cross-strand KRI to the full transcript set.
///
/// Returns (filtered_transcripts, n_killed, n_would_fire).
pub fn cross_strand_kri_filter(
    transcripts: Vec<Transcript>,
    verbose: bool,
) -> (Vec<Transcript>, usize, usize) {
    let dry_run = std::env::var_os("RUSTLE_CROSS_STRAND_DRY_RUN").is_some();
    let min_cov: f64 = std::env::var("RUSTLE_CROSS_STRAND_MIN_COV")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_MIN_LOWINTRON_COV);
    let debug = std::env::var_os("RUSTLE_CROSS_STRAND_DEBUG").is_some();

    let clusters = group_by_overlap(&transcripts);
    let mut dead = vec![false; transcripts.len()];
    let mut n_killed = 0usize;
    let mut n_would_fire = 0usize;

    for cluster in &clusters {
        if cluster.len() < 2 {
            continue;
        }
        // Sort indices by cov desc; ST's predord sorts by cov.
        let mut order: Vec<usize> = cluster.clone();
        order.sort_unstable_by(|&a, &b| {
            transcripts[b]
                .coverage
                .partial_cmp(&transcripts[a].coverage)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Precompute lowintron heuristic per cluster member (cheap).
        let lowintron: Vec<Vec<bool>> = order
            .iter()
            .map(|&i| assumed_lowintron(&transcripts[i], min_cov))
            .collect();

        for ii in 0..order.len() {
            let n1 = order[ii];
            if dead[n1] {
                continue;
            }
            if lowintron[ii].iter().all(|&b| !b) {
                continue; // n1 cov below min_cov — skip
            }
            for jj in (ii + 1)..order.len() {
                let n2 = order[jj];
                if dead[n2] {
                    continue;
                }
                let t1 = &transcripts[n1];
                let t2 = &transcripts[n2];
                if t1.chrom != t2.chrom {
                    continue;
                }
                // Don't kill guide-pinned predictions.
                if t2.ref_transcript_id.is_some() {
                    continue;
                }
                // ST's `overlaps[n1, n2]` gate (rlink.cpp:18957) — only
                // consider pairs that share at least one base of exonic
                // territory. Span-only overlap from cluster grouping is
                // not sufficient.
                if !exons_share_a_base(t1, t2) {
                    continue;
                }
                // SAME-strand pairs are already handled by rustle's
                // per-strand predcluster (which runs Rule A / KRI etc.
                // within each bundle). The whole point of this pass is to
                // catch CROSS-strand kills the per-strand scope misses.
                // Without this gate we step on per-strand predcluster's
                // territory and over-kill same-strand alternative isoforms.
                let allow_same_strand = std::env::var_os(
                    "RUSTLE_CROSS_STRAND_INCLUDE_SAME_STRAND",
                ).is_some();
                if !allow_same_strand && t1.strand == t2.strand {
                    continue;
                }
                let score = retained_intron_score(t1, t2, &lowintron[ii]);
                if score > 0 {
                    n_would_fire += 1;
                    if debug {
                        eprintln!(
                            "[XSKRI] {} n2={}:{}-{}({}) cov={:.2} nex={} \
                             by n1={}:{}-{}({}) cov={:.2} nex={} score={}",
                            if dry_run { "would_kill" } else { "KILL" },
                            t2.chrom,
                            t2.exons.first().map(|e| e.0).unwrap_or(0),
                            t2.exons.last().map(|e| e.1).unwrap_or(0),
                            t2.strand,
                            t2.coverage,
                            t2.exons.len(),
                            t1.chrom,
                            t1.exons.first().map(|e| e.0).unwrap_or(0),
                            t1.exons.last().map(|e| e.1).unwrap_or(0),
                            t1.strand,
                            t1.coverage,
                            t1.exons.len(),
                            score,
                        );
                    }
                    if !dry_run {
                        dead[n2] = true;
                        n_killed += 1;
                    }
                }
            }
        }
    }

    if verbose && (n_killed > 0 || n_would_fire > 0) {
        eprintln!(
            "    cross_strand_kri_filter: would_fire={} killed={} (dry_run={}, min_cov={})",
            n_would_fire, n_killed, dry_run, min_cov
        );
    }

    let filtered: Vec<Transcript> = transcripts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !dead[*i])
        .map(|(_, t)| t)
        .collect();
    (filtered, n_killed, n_would_fire)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(chrom: &str, strand: char, cov: f64, exons: Vec<(u64, u64)>) -> Transcript {
        let mut t = Transcript::default();
        t.chrom = chrom.to_string();
        t.strand = strand;
        t.coverage = cov;
        t.exons = exons;
        t
    }

    #[test]
    fn group_by_overlap_two_strand_neighbors() {
        let txs = vec![
            mk("chr1", '+', 23.0, vec![(100, 200), (300, 400)]),
            mk("chr1", '-', 182.0, vec![(50, 150), (160, 250)]),
            mk("chr2", '+', 5.0, vec![(100, 200)]),
        ];
        let groups = group_by_overlap(&txs);
        assert_eq!(groups.len(), 2);
        let mut g0 = groups[0].clone();
        g0.sort();
        assert_eq!(g0, vec![0, 1]); // chr1 cross-strand pair
        assert_eq!(groups[1], vec![2]); // chr2 alone
    }

    #[test]
    fn assumed_lowintron_threshold() {
        let high = mk("c", '+', 100.0, vec![(0, 10), (20, 30)]);
        let low = mk("c", '+', 1.0, vec![(0, 10), (20, 30)]);
        assert_eq!(assumed_lowintron(&high, 50.0), vec![true]);
        assert_eq!(assumed_lowintron(&low, 50.0), vec![false]);
    }
}
