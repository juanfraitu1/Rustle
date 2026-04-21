//! Junction correction helpers that mirror the original algorithm's inline junction handling.
//! Rustle intentionally does not keep a separate ratio-cluster rewrite pass here,
//! because `` performs the relevant demotions inline in the `higherr` block.

use crate::types::{
    BundleRead, DetHashMap as HashMap, DetHashSet as HashSet, Junction, JunctionStat, JunctionStats,
};

/// junc_ratio_threshold tiers
#[inline]
pub fn junc_ratio_threshold(na: f64) -> f64 {
    if na <= 25.0 {
        1.0
    } else if na <= 100.0 {
        0.50
    } else if na <= 200.0 {
        0.25
    } else {
        0.10
    }
}

/// ok_to_demote gate: current support relative to better support.
#[inline]
pub fn ok_to_demote(cur: &JunctionStat, better: &JunctionStat) -> bool {
    let nb = better.mrcount.max(0.0);
    if nb <= 0.0 {
        return false;
    }
    let na = cur.mrcount.max(0.0);
    (na / nb) <= junc_ratio_threshold(na)
}

/// Mismatch/deletion/insertion near splice anchors (mismatch_anchor-style).
///
/// `anchor_ref_start` follows the original algorithm's `mismatch_anchor(..., refstart, ...)` call site,
/// where `refstart` is the current bundle start.
pub fn mismatch_anchor(read: &BundleRead, junction_support: u64, anchor_ref_start: u64) -> bool {
    if read.exons.len() < 2 {
        return false;
    }
    let seg_len = |idx: usize| -> u64 { read.exons[idx].1.saturating_sub(read.exons[idx].0) };

    // MD pass: mirrors mismatch_anchor behavior, including bundle-start anchoring.
    if let Some(md) = &read.md {
        let bytes = md.as_bytes();
        let mut p = 0usize;
        let mut i = 0usize;
        let mut parsed = anchor_ref_start;
        let mut rdlen = 0u64;

        while p < bytes.len() {
            let c = bytes[p] as char;

            if c.is_ascii_digit() {
                let mut n = 0u64;
                while p < bytes.len() && (bytes[p] as char).is_ascii_digit() {
                    n = n
                        .saturating_mul(10)
                        .saturating_add((bytes[p] - b'0') as u64);
                    p += 1;
                }
                parsed = parsed.saturating_add(n);
                continue;
            }

            if c == '^' {
                p += 1;
                let mut del_len = 0u64;
                while p < bytes.len() && (bytes[p] as char).is_ascii_alphabetic() {
                    del_len = del_len.saturating_add(1);
                    p += 1;
                }

                while i < read.exons.len() && rdlen.saturating_add(seg_len(i)) < parsed {
                    rdlen = rdlen.saturating_add(seg_len(i));
                    i += 1;
                }
                if i == read.exons.len() {
                    break;
                }

                let seg_start = read.exons[i].0;
                let seg_end_excl = read.exons[i].1;
                if (i > 0 && parsed.saturating_sub(seg_start) < junction_support)
                    || (i < read.exons.len() - 1
                        && seg_end_excl.saturating_sub(parsed).saturating_sub(del_len)
                            < junction_support)
                {
                    return true;
                }

                parsed = parsed.saturating_add(del_len);
                continue;
            }

            if c.is_ascii_alphabetic() {
                while i < read.exons.len() && rdlen.saturating_add(seg_len(i)) < parsed {
                    rdlen = rdlen.saturating_add(seg_len(i));
                    i += 1;
                }
                if i == read.exons.len() {
                    break;
                }

                let seg_start = read.exons[i].0;
                let seg_end_excl = read.exons[i].1;
                if (i > 0 && parsed.saturating_sub(seg_start) < junction_support)
                    || (i < read.exons.len() - 1
                        && seg_end_excl.saturating_sub(parsed.saturating_add(1)) < junction_support)
                {
                    return true;
                }

                parsed = parsed.saturating_add(1);
            }

            p += 1;
        }
    }

    // CIGAR insertion pass: mimic behavior (parsed starts at 0 and advances on
    // reference-consuming ops in the original CIGAR scan).
    let mut i = 0usize;
    let mut parsed;
    let mut rdlen = 0u64;
    for ins_abs in &read.insertion_sites {
        parsed = ins_abs.saturating_sub(read.ref_start);
        while i < read.exons.len() && rdlen.saturating_add(seg_len(i)) < parsed {
            rdlen = rdlen.saturating_add(seg_len(i));
            i += 1;
        }
        if i == read.exons.len() {
            break;
        }
        let seg_start = read.exons[i].0;
        let seg_end_excl = read.exons[i].1;
        if (i > 0 && parsed.saturating_sub(seg_start) < junction_support)
            || (i < read.exons.len() - 1
                && seg_end_excl.saturating_sub(parsed.saturating_add(1)) < junction_support)
        {
            return true;
        }
    }

    false
}

/// High-error junction correction: redirect all-bad nearby junctions to a stronger nearby one
/// (subset of higherr branch).
pub fn correct_bundle_junctions_higherr(
    stats: &JunctionStats,
    window: u64,
    tolerance: f64,
    verbose: bool,
) -> (JunctionStats, HashMap<Junction, Junction>) {
    let mut redirect: HashMap<Junction, Junction> = Default::default();
    if stats.is_empty() {
        return (Default::default(), redirect);
    }
    let mut junctions: Vec<Junction> = stats.keys().copied().collect();
    junctions.sort_by_key(|j| (j.donor, j.acceptor));

    // : nm==nreads marks a junction as having ONLY
    // mismatch-bearing reads. For long reads, nm==nreads is ALWAYS true (every
    // long read counts as mismatch in processRead:915). the original implementation only kills
    // such junctions if they're long introns (>100kb) with very low coverage
    // (<10 reads). It does NOT redirect them to nearby stronger junctions.
    //
    // The previous is_bad check (nm >= mrcount) was too aggressive for long reads:
    // it flagged ALL long-read-only junctions as bad, redirecting rare but valid
    // junctions to nearby strong ones. This killed small exons that the original implementation keeps.
    //
    // Fix: require nm > mrcount (strict inequality) OR require the junction to be
    // a long intron with low coverage, matching the original good_junc logic.
    let longintron = 100_000u64;
    let chi_win_error = 10.0f64; // CHI_WIN(100) * ERROR_PERC(0.1)
    let is_bad = |s: &JunctionStat, j: &Junction| {
        if s.nm <= 0.0 || s.mrcount <= 0.0 {
            return false;
        }
        if s.nm < s.mrcount {
            return false; // some reads are mismatch-free → not bad
        }
        // nm >= mrcount: all reads have mismatches. For long reads this is always true.
        // Only flag as bad if it's a long intron with very low coverage (
        let intron_len = j.acceptor.saturating_sub(j.donor);
        intron_len > longintron && s.nreads_good < chi_win_error
    };
    let snap_strict = std::env::var_os("RUSTLE_SNAP_STRICT").is_some();
    let is_reliable = |s: &JunctionStat, j: &Junction| {
        if !(s.mrcount > 0.0 && !is_bad(s, j)) {
            return false;
        }
        if snap_strict && s.nreads_good <= 0.0 {
            return false;
        }
        true
    };

    let mut donor_redirect: HashMap<Junction, Junction> = Default::default();
    let mut acceptor_redirect: HashMap<Junction, Junction> = Default::default();

    fn resolve_redirect(mut j: Junction, redir: &HashMap<Junction, Junction>) -> Junction {
        let mut seen: HashSet<Junction> = Default::default();
        while let Some(nxt) = redir.get(&j).copied() {
            if nxt == j || !seen.insert(j) {
                break;
            }
            j = nxt;
        }
        j
    }

    for j in &junctions {
        let Some(cur) = stats.get(j) else { continue };
        if !is_bad(cur, j) {
            continue;
        }
        let mut best: Option<(Junction, f64)> = None;
        for cand_j in &junctions {
            if cand_j == j || j.donor.abs_diff(cand_j.donor) > window {
                continue;
            }
            let Some(cand_raw) = stats.get(cand_j) else {
                continue;
            };
            if let (Some(a), Some(b)) = (cur.strand, cand_raw.strand) {
                if a != b {
                    continue;
                }
            }
            let cand_j = resolve_redirect(*cand_j, &donor_redirect);
            let Some(cand) = stats.get(&cand_j) else {
                continue;
            };
            if !is_reliable(cand, &cand_j) {
                continue;
            }
            if !(cand.leftsupport > cur.leftsupport * tolerance) {
                continue;
            }
            if !ok_to_demote(cur, cand) {
                continue;
            }
            let score = cand.leftsupport + cand.mrcount;
            match best {
                Some((_, s)) if score <= s => {}
                _ => best = Some((cand_j, score)),
            }
        }
        if let Some((target, _)) = best {
            donor_redirect.insert(*j, target);
        }
    }

    for j in &junctions {
        let Some(cur) = stats.get(j) else { continue };
        if !is_bad(cur, j) {
            continue;
        }
        let mut best: Option<(Junction, f64)> = None;
        for cand_j in &junctions {
            if cand_j == j || j.acceptor.abs_diff(cand_j.acceptor) > window {
                continue;
            }
            let Some(cand_raw) = stats.get(cand_j) else {
                continue;
            };
            if let (Some(a), Some(b)) = (cur.strand, cand_raw.strand) {
                if a != b {
                    continue;
                }
            }
            let cand_j = resolve_redirect(*cand_j, &acceptor_redirect);
            let Some(cand) = stats.get(&cand_j) else {
                continue;
            };
            if !is_reliable(cand, &cand_j) {
                continue;
            }
            if !(cand.rightsupport > cur.rightsupport * tolerance) {
                continue;
            }
            if !ok_to_demote(cur, cand) {
                continue;
            }
            let score = cand.rightsupport + cand.mrcount;
            match best {
                Some((_, s)) if score <= s => {}
                _ => best = Some((cand_j, score)),
            }
        }
        if let Some((target, _)) = best {
            acceptor_redirect.insert(*j, target);
        }
    }

    for j in &junctions {
        let Some(cur) = stats.get(j) else { continue };
        if !is_bad(cur, j) {
            continue;
        }
        let left = donor_redirect.get(j).copied();
        let right = acceptor_redirect.get(j).copied();
        let target = match (left, right) {
            (Some(l), Some(r)) if l == r => Some(l),
            (Some(l), Some(r)) => {
                let synth = Junction::new(l.donor, r.acceptor);
                if let Some(s) = stats.get(&synth) {
                    if is_reliable(s, &synth) {
                        Some(synth)
                    } else {
                        Some(l)
                    }
                } else {
                    Some(l)
                }
            }
            (Some(l), None) => Some(l),
            (None, Some(r)) => Some(r),
            _ => None,
        };
        if let Some(t) = target {
            if t != *j {
                redirect.insert(*j, t);
            }
        }
    }

    let mut new_stats: JunctionStats = Default::default();
    for (j, stat) in stats.iter() {
        let canonical = redirect.get(j).copied().unwrap_or(*j);
        let entry = new_stats
            .entry(canonical)
            .or_insert_with(JunctionStat::default);
        entry.mrcount += stat.mrcount;
        entry.nreads_good += stat.nreads_good;
        entry.rcount += stat.rcount;
        entry.nm += stat.nm;
        entry.mm += stat.mm;
        entry.leftsupport += stat.leftsupport;
        entry.rightsupport += stat.rightsupport;
        if entry.strand.is_none() {
            entry.strand = stat.strand;
        }
    }

    if std::env::var_os("RUSTLE_GOODJUNC_DEBUG").is_some() {
        for (from, to) in &redirect {
            eprintln!(
                "HE_REDIRECT {}-{} → {}-{} donor_dist={} acceptor_dist={}",
                from.donor,
                from.acceptor,
                to.donor,
                to.acceptor,
                from.donor.abs_diff(to.donor),
                from.acceptor.abs_diff(to.acceptor)
            );
        }
    }
    if verbose && !redirect.is_empty() {
        eprintln!(
            "    higherr correction redirected {} junction(s)",
            redirect.len()
        );
    }
    (new_stats, redirect)
}

/// Demote alt-donor and alt-acceptor junctions that have much less read
/// support than a nearby canonical junction sharing the other endpoint.
/// Returns redirect map: alt_junction → canonical_junction.
///
/// Reduces graph segmentation at alt-splice hot spots (STRG.309 class):
/// after reads get their alt-donor boundaries snapped to canonical, the
/// resulting graph has fewer distinguishing segments, so max_flow explores
/// fewer path combinations and emits fewer noise variants.
///
/// Targets ONLY intra-acceptor (or intra-donor) alt-splicing — does not
/// touch junctions whose donor AND acceptor both differ from canonicals.
/// Those are true alt-splice events that should be preserved.
///
/// Env vars (all opt-in):
///   RUSTLE_ALT_JUNC_SNAP=1       enable
///   RUSTLE_ALT_JUNC_SNAP_WINDOW  bp tolerance (default 50)
///   RUSTLE_ALT_JUNC_SNAP_RATIO   alt must have at least this fraction of
///                                 canonical reads to survive (default 0.5)
///   RUSTLE_ALT_JUNC_SNAP_DEBUG=1 trace each demotion
pub fn demote_alt_junctions_to_canonical(
    stats: &JunctionStats,
) -> HashMap<Junction, Junction> {
    let mut redirect: HashMap<Junction, Junction> = Default::default();
    if std::env::var_os("RUSTLE_ALT_JUNC_SNAP").is_none() {
        return redirect;
    }
    let window: u64 = std::env::var("RUSTLE_ALT_JUNC_SNAP_WINDOW")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(50);
    let min_ratio: f64 = std::env::var("RUSTLE_ALT_JUNC_SNAP_RATIO")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0.5);
    let debug = std::env::var_os("RUSTLE_ALT_JUNC_SNAP_DEBUG").is_some();

    // Skip junctions marked as guide (preserve explicit annotation)
    let reads_of = |j: &Junction| -> f64 {
        stats.get(j).map(|s| {
            if s.guide_match { return f64::MAX; } // guide junctions never demoted
            s.mrcount.max(s.nreads_good)
        }).unwrap_or(0.0)
    };

    let juncs: Vec<&Junction> = stats.keys().collect();

    // 1. Alt-donor demotion: for each acceptor, find canonical donor, demote
    //    alt-donors within WINDOW with < MIN_RATIO * canonical reads.
    use std::collections::HashMap as StdHashMap;
    let mut by_acceptor: StdHashMap<u64, Vec<&Junction>> = StdHashMap::new();
    for &j in &juncs {
        by_acceptor.entry(j.acceptor).or_default().push(j);
    }
    for (_acc, sibs) in &by_acceptor {
        if sibs.len() < 2 { continue; }
        let canon = sibs.iter().copied().max_by(|a, b| {
            reads_of(a).partial_cmp(&reads_of(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(canon) = canon else { continue; };
        let canon_r = reads_of(canon);
        if canon_r <= 0.0 || canon_r == f64::MAX { continue; }
        for &alt in sibs {
            if alt.donor == canon.donor { continue; }
            // Skip guide alt-junctions
            if stats.get(alt).map(|s| s.guide_match).unwrap_or(false) { continue; }
            let alt_r = reads_of(alt);
            if alt_r == f64::MAX { continue; }
            let donor_dist = alt.donor.abs_diff(canon.donor);
            if donor_dist <= window && alt_r < min_ratio * canon_r {
                redirect.insert(*alt, *canon);
                if debug {
                    eprintln!(
                        "ALT_JUNC_SNAP demote_at_acceptor={} alt={}-{} reads={:.0} → canon={}-{} reads={:.0} dist={}bp",
                        alt.acceptor, alt.donor, alt.acceptor, alt_r,
                        canon.donor, canon.acceptor, canon_r, donor_dist
                    );
                }
            }
        }
    }

    // 2. Alt-acceptor demotion: for each donor, find canonical acceptor,
    //    demote alt-acceptors within WINDOW with < MIN_RATIO * canonical reads.
    let mut by_donor: StdHashMap<u64, Vec<&Junction>> = StdHashMap::new();
    for &j in &juncs {
        by_donor.entry(j.donor).or_default().push(j);
    }
    for (_don, sibs) in &by_donor {
        if sibs.len() < 2 { continue; }
        let canon = sibs.iter().copied().max_by(|a, b| {
            reads_of(a).partial_cmp(&reads_of(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(canon) = canon else { continue; };
        let canon_r = reads_of(canon);
        if canon_r <= 0.0 || canon_r == f64::MAX { continue; }
        for &alt in sibs {
            if alt.acceptor == canon.acceptor { continue; }
            if stats.get(alt).map(|s| s.guide_match).unwrap_or(false) { continue; }
            let alt_r = reads_of(alt);
            if alt_r == f64::MAX { continue; }
            let acc_dist = alt.acceptor.abs_diff(canon.acceptor);
            if acc_dist <= window && alt_r < min_ratio * canon_r {
                // Don't overwrite existing redirect (alt-donor takes priority)
                redirect.entry(*alt).or_insert(*canon);
                if debug {
                    eprintln!(
                        "ALT_JUNC_SNAP demote_at_donor={} alt={}-{} reads={:.0} → canon={}-{} reads={:.0} dist={}bp",
                        alt.donor, alt.donor, alt.acceptor, alt_r,
                        canon.donor, canon.acceptor, canon_r, acc_dist
                    );
                }
            }
        }
    }

    if debug && !redirect.is_empty() {
        eprintln!(
            "ALT_JUNC_SNAP total={} window={}bp ratio={:.2}",
            redirect.len(), window, min_ratio
        );
    }
    redirect
}

/// Guide-based coordinate snapping for high-mismatch junctions
/// For junctions where nm >= nreads (all reads are high-mismatch), snaps donor/acceptor
/// to the closest guide junction within `sserror` tolerance.
/// Returns (updated JunctionStats, correction map from old Junction to new Junction).
pub fn snap_junctions_to_guides(
    stats: &JunctionStats,
    guide_junctions: &[Junction],
    sserror: u64,
    longreads: bool,
) -> (JunctionStats, HashMap<Junction, Junction>) {
    let mut snap_map: HashMap<Junction, Junction> = Default::default();
    if stats.is_empty() || guide_junctions.is_empty() {
        return (stats.clone(), snap_map);
    }

    let gjd = std::env::var_os("RUSTLE_GOODJUNC_DEBUG").is_some();

    // Build sorted guide junction lists by donor and by acceptor
    let mut guides_by_donor: Vec<Junction> = guide_junctions.to_vec();
    guides_by_donor.sort_by_key(|j| j.donor);
    guides_by_donor.dedup();
    let mut guides_by_acceptor: Vec<Junction> = guide_junctions.to_vec();
    guides_by_acceptor.sort_by_key(|j| j.acceptor);
    guides_by_acceptor.dedup();

    // Sort junctions by donor (junction list is sorted by start)
    let mut junctions: Vec<Junction> = stats.keys().copied().collect();
    junctions.sort_by_key(|j| (j.donor, j.acceptor));

    // Sort the same junctions by acceptor (ejunction sorted by end)
    let mut junctions_by_acc: Vec<Junction> = junctions.clone();
    junctions_by_acc.sort_by_key(|j| (j.acceptor, j.donor));

    // Build index maps: junction → its index in each sorted order
    let mut donor_idx_map: HashMap<Junction, usize> = Default::default();
    for (i, j) in junctions.iter().enumerate() {
        donor_idx_map.insert(*j, i);
    }
    let mut acc_idx_map: HashMap<Junction, usize> = Default::default();
    for (i, j) in junctions_by_acc.iter().enumerate() {
        acc_idx_map.insert(*j, i);
    }

    // Track which junctions get snapped (donor or acceptor modified)
    // We modify in-place conceptually, tracking old→new mapping
    let mut new_donors: HashMap<Junction, u64> = Default::default();
    let mut new_acceptors: HashMap<Junction, u64> = Default::default();
    let _new_strands: HashMap<Junction, Option<i8>> = Default::default();

    // Pass 1: snap donors (lines 16710-16750)
    // Use a sliding window pointer
    let mut s = 0usize;
    for jn in &junctions {
        let st = match stats.get(jn) {
            Some(s) => s,
            None => continue,
        };
        if st.guide_match {
            continue;
        }
        if !(st.nm > 0.0 && st.nm.round() >= st.mrcount.round()) {
            continue;
        }

        // Advance s past guides whose donor + sserror < jn.donor
        while s < guides_by_donor.len() && guides_by_donor[s].donor + sserror < jn.donor {
            s += 1;
        }

        let mut k = s;
        let mut best_idx: Option<usize> = None;
        let mut best_dist = sserror + 1;

        while k < guides_by_donor.len() && guides_by_donor[k].donor <= jn.donor + sserror {
            let gj = &guides_by_donor[k];
            // strand check: (!junction[i]->strand || gjunc[k]->strand==junction[i]->strand || longreads)
            // In longreads mode, always accept
            if longreads || st.strand.is_none() || st.strand == Some(0) {
                // Accept any guide
            } else {
                // Would need guide strand info, but our guide junctions don't carry strand.
                // In longreads mode (our case), this branch is skipped.
            }

            if gj.donor == jn.donor {
                // Perfect match on donor — no need to snap
                best_idx = None;
                break;
            }

            let d = gj.donor.abs_diff(jn.donor);
            if d < best_dist {
                best_dist = d;
                best_idx = Some(k);
            }
            k += 1;
        }

        if let Some(bi) = best_idx {
            let new_donor = guides_by_donor[bi].donor;
            if gjd {
                eprintln!(
                    "GUIDE_SNAP_DONOR {}-{} donor {} → {} dist={}",
                    jn.donor, jn.acceptor, jn.donor, new_donor, best_dist
                );
            }
            new_donors.insert(*jn, new_donor);
        }
    }

    // Pass 2: snap acceptors (lines 16751-16789)
    let mut e = 0usize;
    for jn in &junctions_by_acc {
        let st = match stats.get(jn) {
            Some(s) => s,
            None => continue,
        };
        if st.guide_match {
            continue;
        }
        if !(st.nm > 0.0 && st.nm.round() >= st.mrcount.round()) {
            continue;
        }

        // Advance e past guides whose acceptor + sserror < jn.acceptor
        while e < guides_by_acceptor.len() && guides_by_acceptor[e].acceptor + sserror < jn.acceptor
        {
            e += 1;
        }

        let mut k = e;
        let mut best_idx: Option<usize> = None;
        let mut best_dist = sserror + 1;

        while k < guides_by_acceptor.len()
            && guides_by_acceptor[k].acceptor <= jn.acceptor + sserror
        {
            let gj = &guides_by_acceptor[k];

            if gj.acceptor == jn.acceptor {
                // Perfect match on acceptor — no need to snap
                best_idx = None;
                break;
            }

            let d = gj.acceptor.abs_diff(jn.acceptor);
            if d < best_dist {
                best_dist = d;
                best_idx = Some(k);
            }
            k += 1;
        }

        if let Some(bi) = best_idx {
            let new_acc = guides_by_acceptor[bi].acceptor;
            if gjd {
                eprintln!(
                    "GUIDE_SNAP_ACCEPTOR {}-{} acceptor {} → {} dist={}",
                    jn.donor, jn.acceptor, jn.acceptor, new_acc, best_dist
                );
            }
            new_acceptors.insert(*jn, new_acc);
        }
    }

    // Build the final snapped junctions and correction map
    let mut new_stats: JunctionStats = Default::default();
    for (jn, st) in stats.iter() {
        let snapped_donor = new_donors.get(jn).copied().unwrap_or(jn.donor);
        let snapped_acc = new_acceptors.get(jn).copied().unwrap_or(jn.acceptor);
        let new_jn = Junction::new(snapped_donor, snapped_acc);

        if new_jn != *jn {
            snap_map.insert(*jn, new_jn);
        }

        // Merge into new stats (deduplication: junctions snapped to same coords merge)
        let entry = new_stats
            .entry(new_jn)
            .or_insert_with(JunctionStat::default);
        if entry.mrcount == 0.0 {
            // First entry for this junction — copy all fields
            *entry = st.clone();
        } else {
            // Merge: sum counts
            entry.mrcount += st.mrcount;
            entry.nm += st.nm;
            entry.mm += st.mm;
            entry.nreads_good += st.nreads_good;
            entry.leftsupport += st.leftsupport;
            entry.rightsupport += st.rightsupport;
        }

        // Check if snapped junction now matches a guide exactly
        if new_jn != *jn {
            let exact_match = guide_junctions
                .iter()
                .any(|gj| gj.donor == new_jn.donor && gj.acceptor == new_jn.acceptor);
            if exact_match {
                if let Some(e) = new_stats.get_mut(&new_jn) {
                    e.guide_match = true;
                }
            }
        }
    }

    // Kill junctions that were merged into another (sets strand=0, nreads=0)
    // These are junctions whose old key still exists but got merged into a different new key.
    // The snap_map already tracks old→new, so reads can be remapped.
    // We also need to handle the case where the original junction key is NOT in new_stats
    // (because it was remapped) — that's fine, those entries merged into the new key.

    if gjd && !snap_map.is_empty() {
        eprintln!("GUIDE_SNAP: {} junctions snapped", snap_map.len());
    }

    (new_stats, snap_map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::JunctionStat;

    #[test]
    fn test_higherr_merges_weak_nearby() {
        let mut stats: JunctionStats = Default::default();
        // Strong junction at (100, 200)
        stats.insert(
            Junction::new(100, 200),
            JunctionStat {
                mrcount: 10.0,
                leftsupport: 10.0,
                rightsupport: 10.0,
                rcount: 10,
                ..Default::default()
            },
        );
        // Weak nearby at (105, 205) - within 25bp
        stats.insert(
            Junction::new(105, 205),
            JunctionStat {
                mrcount: 1.0,
                nm: 1.0,
                leftsupport: 1.0,
                rightsupport: 1.0,
                rcount: 1,
                ..Default::default()
            },
        );
        let (new_stats, map) = correct_bundle_junctions_higherr(&stats, 25, 0.9, false);
        assert!(map.contains_key(&Junction::new(105, 205)));
        assert_eq!(
            map.get(&Junction::new(105, 205)),
            Some(&Junction::new(100, 200))
        );
        assert_eq!(new_stats.len(), 1);
        assert!((new_stats.get(&Junction::new(100, 200)).unwrap().mrcount - 11.0).abs() < 1e-6);
    }
}
