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

/// Total number of bases shared between exons of `a` and `b`.
fn exon_overlap_bp(a: &Transcript, b: &Transcript) -> u64 {
    let mut i = 0;
    let mut j = 0;
    let mut total: u64 = 0;
    while i < a.exons.len() && j < b.exons.len() {
        let (as_, ae) = a.exons[i];
        let (bs, be) = b.exons[j];
        let s = as_.max(bs);
        let e = ae.min(be);
        if s < e {
            total += e - s;
        }
        if ae <= be {
            i += 1;
        } else {
            j += 1;
        }
    }
    total
}

fn transcript_tlen(t: &Transcript) -> u64 {
    t.exons.iter().map(|&(s, e)| e.saturating_sub(s)).sum()
}

/// Mirror ST's `update_overlap` per-engulfing-clause threshold
/// (rlink.cpp:17828-17855). For each of four engulfing configurations,
/// reject the overlap when the engulfed-region length is shorter than
/// `error_perc` of the engulfing transcript's tlen.
///
/// The four clauses (using ST's pred-pair p, pi):
///   1. p's first exon engulfs pi: `p.exons[0].end >= pi.end`
///      → len = pi.end - p.start + 1, threshold = error_perc * p.tlen
///   2. pi's first exon engulfs p: `pi.exons[0].end >= p.end`
///      → len = p.end - pi.start + 1, threshold = error_perc * pi.tlen
///   3. p's last exon contains pi.start: `p.exons.last.start <= pi.start`
///      → len = p.end - pi.start + 1, threshold = error_perc * p.tlen
///   4. pi's last exon contains p.start: `pi.exons.last.start <= p.start`
///      → len = pi.end - p.start + 1, threshold = error_perc * pi.tlen
///
/// We treat the pair (a, b) symmetrically: try each clause with both
/// orderings (a as p, b as pi). The overlap is significant only if NO
/// matching clause rejects.
///
/// Crucially (vs the simpler min/max heuristics): clause 2 uses pi.tlen
/// — the LARGER tlen when pi engulfs p — which correctly rejects the
/// STRG.136 case (pi=STRG.136.4 first exon engulfs p=RSTL.366.1's end;
/// 119 bp < 10% × 3932 = 393), while keeping the RSTL.135 case (overlap
/// 1091 bp > 10% × 2866 = 286 from pi=RSTL.135.14's first exon).
fn exons_significantly_overlap(
    a: &Transcript,
    b: &Transcript,
    error_perc: f64,
) -> bool {
    if exon_overlap_bp(a, b) == 0 {
        return false;
    }
    !any_clause_rejects(a, b, error_perc) && !any_clause_rejects(b, a, error_perc)
}

/// Run all four ST update_overlap clauses with `p=a, pi=b`. Returns true
/// if any clause's len < threshold (i.e. ST would set overlaps=false).
/// Caller flips arguments to test the symmetric case.
fn any_clause_rejects(p: &Transcript, pi: &Transcript, error_perc: f64) -> bool {
    let p_exons = &p.exons;
    let pi_exons = &pi.exons;
    if p_exons.is_empty() || pi_exons.is_empty() {
        return false;
    }
    let p_start = p_exons[0].0;
    let p_end = p_exons.last().unwrap().1;
    let pi_start = pi_exons[0].0;
    let pi_end = pi_exons.last().unwrap().1;
    let p_tlen = transcript_tlen(p) as f64;
    let pi_tlen = transcript_tlen(pi) as f64;
    let p_first_end = p_exons[0].1;
    let p_last_start = p_exons.last().unwrap().0;

    // Clause 1: p's first exon engulfs pi entirely (`p.exons[0].end >= pi.end`).
    if p_first_end >= pi_end {
        let len = pi_end.saturating_sub(p_start) + 1;
        if (len as f64) < error_perc * p_tlen {
            return true;
        }
    }
    // Clause 3: p's last exon contains pi's start
    // (`p.exons.last.start <= pi.start`).
    if p_last_start <= pi_start {
        let len = p_end.saturating_sub(pi_start) + 1;
        if (len as f64) < error_perc * p_tlen {
            return true;
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
                // ST's `overlaps[n1, n2]` gate (rlink.cpp:18957) plus the
                // `update_overlap` size-threshold rejection
                // (rlink.cpp:17828-17855). Skips antisense-overlapping-
                // termini pairs whose shared region is < ERROR_PERC ×
                // min(tlen).
                if !exons_significantly_overlap(t1, t2, 0.1) {
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

    #[test]
    fn significant_overlap_rejects_strg_136_pattern() {
        // STRG.136 case: p (4 exons, tlen 974) ends at 23081962; pi
        // (19 exons, tlen 3932) first exon = 23081844-23083200 engulfs p's end.
        // Clause 2 (pi.exons[0].end >= p.end): len = 23081962 - 23081844 + 1 = 119.
        // 119 < 0.1 × 3932 = 393 → reject.
        let p = mk(
            "c",
            '-',
            66.0,
            vec![
                (23080202, 23080591),
                (23080869, 23081007),
                (23081364, 23081434),
                (23081589, 23081962),
            ],
        );
        let pi = mk("c", '+', 3.3, vec![(23081844, 23083200), (23083988, 23085000)]);
        // Use simplified pi tlen to keep the test focused: 1357 + 1013 = 2370.
        // Clause 2: 119 bp / 0.1 × 2370 = 119 / 237 → REJECT.
        assert!(!exons_significantly_overlap(&p, &pi, 0.1));
    }

    #[test]
    fn significant_overlap_keeps_rstl_135_pattern() {
        // RSTL.135 case: p=- strand antisense (cov 174, 11 exons, tlen 2180),
        // pi=+ strand RSTL.135.14 first exon engulfs p's end with 1091 bp.
        // Clause 2: len = 1091, threshold = 0.1 × pi.tlen (2866) = 286.6.
        // 1091 > 286.6 → keep.
        let p = mk(
            "c",
            '-',
            174.0,
            vec![
                (44646789, 44647096),
                (44648168, 44648350),
                (44649248, 44649329),
                (44651129, 44651310),
                (44654733, 44654922),
                (44657248, 44657438),
                (44661730, 44662177),
                (44664021, 44664141),
                (44664287, 44664387),
                (44666408, 44666510),
                (44670562, 44670832),
            ],
        );
        // RSTL.135.14: 11 exons spanning to 44715975, tlen ≈ 2866.
        let pi = mk(
            "c",
            '+',
            23.8,
            vec![
                (44669742, 44670957),
                (44672486, 44672587),
                (44678780, 44678929),
                (44681195, 44681251),
                (44682100, 44682172),
                (44690899, 44690945),
                (44704565, 44704647),
                (44704745, 44704843),
                (44705182, 44705308),
                (44712215, 44712268),
                (44715107, 44715975),
            ],
        );
        assert!(exons_significantly_overlap(&p, &pi, 0.1));
    }

    #[test]
    fn significant_overlap_keeps_real_inclusion() {
        // n2 mostly inside n1 — clear inclusion, should be kept.
        let a = mk("c", '+', 100.0, vec![(0, 1000), (1500, 2500)]);
        let b = mk("c", '+', 10.0, vec![(100, 900)]);
        assert_eq!(exon_overlap_bp(&a, &b), 800);
        assert!(exons_significantly_overlap(&a, &b, 0.1));
    }

}
