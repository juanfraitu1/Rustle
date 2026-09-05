//! Gene-conversion mosaic-read detection (audit theme: "VG finds unusual exon combinations").
//!
//! A paralog family has near-identical copies that differ at a set of *diagnostic sites*.
//! For a read, at each diagnostic site it covers we know which copies' expected base it
//! matches. A GENE-CONVERSION recombinant read's per-site copy pattern SWITCHES from one
//! copy to another at a contiguous breakpoint (copy A for a run of sites, then copy B) —
//! an "unusual combination" directly observed in ONE read, not enumerated. This module
//! detects that switch and rejects the look-alikes: sequencing-error flips, non-identifiable
//! reads, low-power reads, and (at the family layer) one-off chimeras.
//!
//! Design = synthesis of an independent design panel (statistical / algorithmic / biological
//! lenses). The per-read core `detect_mosaic` is a PURE function over the per-site match
//! structure; `aggregate_family` confirms a conversion only when the breakpoint RECURS across
//! independent molecules (the conversion-vs-chimera discriminator). Default-OFF in the
//! pipeline (RUSTLE_VG_MOSAIC_ON); additive — the existing per-copy EM scoring is untouched.
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

/// Per-site observation handed to the detector: which copies the read matches at one
/// diagnostic site. `match_bits[c]` = (read base == copy c's expected base). Sites are in
/// genomic order.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteObs {
    pub ref_pos: u64,
    pub match_bits: Vec<bool>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MosaicStatus {
    LowPower,         // too few decisive sites to call anything
    NonIdentifiable,  // copies tie at too many sites — no per-site signal; abstain
    SingleCopy,       // all decisive sites point at one copy
    NoSwitch,         // a switch was scored but failed the gates (error scatter / weak)
    Mosaic,           // a confident contiguous copy switch
}

/// Per-read verdict. Abstain/SingleCopy/NoSwitch are all `is_mosaic()==false`; the EM's
/// existing per-copy scoring is independent of this (additive metadata).
#[derive(Debug, Clone, PartialEq)]
pub struct MosaicCall {
    pub status: MosaicStatus,
    pub n_sites: usize,
    pub n_decisive: usize,
    pub copy_a: Option<usize>,            // 5'-proximal tract copy
    pub copy_b: Option<usize>,            // 3'-proximal tract copy
    pub breakpoint_ref: Option<(u64, u64)>, // (last A-site, first B-site) — honest bracket
    pub tract_a_sites: usize,
    pub tract_b_sites: usize,
    pub left_purity: f64,
    pub right_purity: f64,
    pub margin: i32,                      // decisive sites the split explains beyond best single
    pub lr_switch: f64,                   // 2*(two-segment loglik − best single-copy loglik)
    pub threshold_used: f64,
    pub score: f64,                       // lr_switch − threshold_used (>0 iff Mosaic)
}

impl MosaicCall {
    pub fn is_mosaic(&self) -> bool {
        self.status == MosaicStatus::Mosaic
    }
}

#[derive(Debug, Clone)]
pub struct MosaicParams {
    pub min_decisive_sites: usize, // D below this → LowPower
    pub min_tract_sites: usize,    // each flank must have ≥ this many agreeing decisive sites
    pub min_improvement: i32,      // split must explain ≥ this many more decisive sites
    pub min_seg_purity: f64,       // per-tract agreement fraction on both sides
    pub max_ambig_frac: f64,       // > this fraction Ambig → NonIdentifiable
    pub alpha_target: f64,         // target per-read FPR (Bonferroni-corrected internally)
    pub bic_penalty_coeff: f64,    // coeff on ln(S) in the model-complexity penalty
    pub eps_floor: f64,            // lower clamp on per-read error rate
    pub eps_cap: f64,              // upper clamp (and fail-safe when read.de is absent)
    // family aggregation
    pub family_min_supporting_reads: usize,
    pub breakpoint_tol: u64,
    pub max_breakpoint_dispersion: u64,
}

impl Default for MosaicParams {
    fn default() -> Self {
        MosaicParams {
            min_decisive_sites: 6,
            min_tract_sites: 3,
            min_improvement: 3,
            min_seg_purity: 0.85,
            max_ambig_frac: 0.5,
            alpha_target: 0.01,
            bic_penalty_coeff: 2.0,
            eps_floor: 0.005,
            eps_cap: 0.05,
            family_min_supporting_reads: 3,
            breakpoint_tol: 50,
            max_breakpoint_dispersion: 50,
        }
    }
}

impl MosaicParams {
    /// Build from RUSTLE_VG_MOSAIC_* env overrides, falling back to defaults.
    pub fn from_env() -> Self {
        let mut p = MosaicParams::default();
        let getf = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
        let getu = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
        if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_DECISIVE") { p.min_decisive_sites = v; }
        if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_TRACT") { p.min_tract_sites = v; }
        if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_IMPROVEMENT") { p.min_improvement = v as i32; }
        if let Some(v) = getf("RUSTLE_VG_MOSAIC_MIN_PURITY") { p.min_seg_purity = v; }
        if let Some(v) = getf("RUSTLE_VG_MOSAIC_ALPHA") { p.alpha_target = v; }
        if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_READS") { p.family_min_supporting_reads = v; }
        p
    }

    /// Per-read error rate clamped to [eps_floor, eps_cap]; `None` → eps_cap (fail-safe).
    pub fn eps_for(&self, de: Option<f32>) -> f64 {
        de.map(|d| (d as f64).clamp(self.eps_floor, self.eps_cap)).unwrap_or(self.eps_cap)
    }
}

/// Pure per-read detector. `obs` = ordered per-site match structure; `eps` = per-site error
/// rate (clamped by the caller). Deterministic, no I/O. With decisive-site agreements the
/// likelihood-ratio reduces to `2·margin·ln((1−eps)/eps)`, and the χ²(2df) null quantile is
/// exactly `−2·ln(α)` — so the threshold is closed-form (the spec's bootstrap_reps=0 path).
pub fn detect_mosaic(obs: &[SiteObs], n_copies: usize, eps: f64, p: &MosaicParams) -> MosaicCall {
    let s = obs.len();
    let mk = |status: MosaicStatus, d: usize| MosaicCall {
        status, n_sites: s, n_decisive: d,
        copy_a: None, copy_b: None, breakpoint_ref: None,
        tract_a_sites: 0, tract_b_sites: 0, left_purity: 0.0, right_purity: 0.0,
        margin: 0, lr_switch: 0.0, threshold_used: 0.0, score: 0.0,
    };

    // Tokenize: a DECISIVE site has exactly one matching copy; ≥2 = Ambig (non-identifiable
    // here); 0 = Novel (read base matches no modeled copy) — a wildcard, ignored.
    let mut decisive_copy: Vec<usize> = Vec::with_capacity(s);
    let mut decisive_pos: Vec<u64> = Vec::with_capacity(s);
    let mut n_ambig = 0usize;
    for o in obs {
        let mut matched = usize::MAX;
        let mut n_match = 0usize;
        for c in 0..n_copies {
            if o.match_bits.get(c).copied().unwrap_or(false) {
                n_match += 1;
                matched = c;
            }
        }
        if n_match == 1 {
            decisive_copy.push(matched);
            decisive_pos.push(o.ref_pos);
        } else if n_match >= 2 {
            n_ambig += 1;
        }
    }
    let d = decisive_copy.len();

    // Gates (abstain before any switch test). NonIdentifiable takes precedence over
    // LowPower: a read swamped by ties has no per-site signal even if it covers many sites
    // (and an all-ambiguous read has d=0, which would otherwise read as LowPower).
    if s > 0 && (n_ambig as f64) > p.max_ambig_frac * s as f64 {
        return mk(MosaicStatus::NonIdentifiable, d);
    }
    if d < p.min_decisive_sites {
        return mk(MosaicStatus::LowPower, d);
    }
    // Distinct copies among decisive sites, ASCENDING (a dense Vec, not a BTreeSet: avoids a
    // per-read heap allocation + tree walk in the hot scan below; ascending order preserves the
    // old tie-break — `max_by_key` keeps the highest copy id among ties).
    let mut distinct: Vec<usize> = decisive_copy.clone();
    distinct.sort_unstable();
    distinct.dedup();
    if distinct.len() <= 1 {
        return mk(MosaicStatus::SingleCopy, d);
    }

    // Prefix counts make a range agreement count O(1): `prefix[c][i]` = #copy-c decisive sites
    // in `decisive_copy[0..i]`. Replaces the per-(k, copy) linear re-scan
    // (O(d²·copies) → O(d·copies)); one allocation reused across all k.
    let ncap = distinct.last().copied().unwrap_or(0) + 1;
    let mut prefix: Vec<Vec<u32>> = vec![vec![0u32; d + 1]; ncap];
    for i in 0..d {
        for c in 0..ncap {
            prefix[c][i + 1] = prefix[c][i];
        }
        prefix[decisive_copy[i]][i + 1] += 1;
    }
    let agree_for = |c: usize, lo: usize, hi: usize| -> usize {
        (prefix[c][hi] - prefix[c][lo]) as usize
    };
    let best_single = distinct.iter().map(|&c| agree_for(c, 0, d)).max().unwrap_or(0);

    // Exhaustive single-changepoint scan over decisive-site split indices, maximizing the
    // number of decisive sites explained by (copy a left of k, copy b right of k), a≠b.
    let mut best_split = best_single;
    let mut best_kab: Option<(usize, usize, usize)> = None;
    for k in 1..d {
        let a = *distinct.iter().max_by_key(|&&c| agree_for(c, 0, k)).unwrap();
        let b = *distinct.iter().max_by_key(|&&c| agree_for(c, k, d)).unwrap();
        if a == b {
            continue;
        }
        let explained = agree_for(a, 0, k) + agree_for(b, k, d);
        // First k to reach a new best (or to match the single-copy best) wins the tie.
        if explained > best_split || (explained == best_split && best_kab.is_none()) {
            best_split = explained;
            best_kab = Some((k, a, b));
        }
    }

    let (k, a, b) = match best_kab {
        Some(x) => x,
        None => return mk(MosaicStatus::SingleCopy, d), // no a≠b split improves on single copy
    };

    let tract_a_sites = agree_for(a, 0, k);
    let tract_b_sites = agree_for(b, k, d);
    let left_purity = tract_a_sites as f64 / k as f64;
    let right_purity = tract_b_sites as f64 / (d - k) as f64;
    let margin = best_split as i32 - best_single as i32;

    // Likelihood-ratio over decisive sites: 2·margin·ln((1−eps)/eps).
    let lr = 2.0 * margin as f64 * ((1.0 - eps) / eps).ln();
    // Threshold = BIC complexity penalty + χ²(2df) quantile at the Bonferroni-corrected α.
    let n_pairs = (n_copies * n_copies.saturating_sub(1) / 2).max(1) as f64;
    let alpha_corr = (p.alpha_target / (s as f64 * n_pairs)).max(1e-12);
    let threshold = p.bic_penalty_coeff * (s as f64).ln() + (-2.0 * alpha_corr.ln());

    let mut call = mk(MosaicStatus::NoSwitch, d);
    call.copy_a = Some(a);
    call.copy_b = Some(b);
    call.breakpoint_ref = Some((decisive_pos[k - 1], decisive_pos[k]));
    call.tract_a_sites = tract_a_sites;
    call.tract_b_sites = tract_b_sites;
    call.left_purity = left_purity;
    call.right_purity = right_purity;
    call.margin = margin;
    call.lr_switch = lr;
    call.threshold_used = threshold;
    call.score = lr - threshold;

    // Dual gate: calibrated LR AND hard integer/run-length backstops (the backstops hold
    // even when the i.i.d. error null is violated by correlated/homopolymer error bursts).
    let pass = a != b
        && tract_a_sites >= p.min_tract_sites
        && tract_b_sites >= p.min_tract_sites
        && margin >= p.min_improvement
        && left_purity >= p.min_seg_purity
        && right_purity >= p.min_seg_purity
        && lr > threshold;
    if pass {
        call.status = MosaicStatus::Mosaic;
    }
    call
}

/// A family-level confirmed (or suspected) gene-conversion event.
#[derive(Debug, Clone, PartialEq)]
pub struct ConversionEvent {
    pub copy_a: usize,
    pub copy_b: usize,
    pub chrom: String,              // chromosome the breakpoint coordinates live on (for the microhomology check)
    pub breakpoint_ref: (u64, u64), // consensus bracket (min last-A, max first-B)
    pub n_supporting_reads: usize,
    pub breakpoint_dispersion: u64, // spread of per-read breakpoint midpoints
    pub confirmed: bool,            // false = ChimeraSuspect
}

/// Family aggregation: a genuine conversion RECURS at a fixed breakpoint across independent
/// molecules; a one-off chimera does not. Inputs are per-read Mosaic calls (caller dedupes
/// to distinct molecules first). Clusters by oriented (a→b) pair and breakpoint midpoint
/// within `breakpoint_tol`; confirms a cluster with ≥ `family_min_supporting_reads` molecules
/// and tight dispersion.
pub fn aggregate_family(calls: &[MosaicCall], chroms: &[&str], p: &MosaicParams) -> Vec<ConversionEvent> {
    // `chroms[i]` is the chromosome `calls[i]`'s breakpoint coordinates live on (parallel array). The
    // chrom is part of the cluster KEY: breakpoints at coincidentally-similar positions on DIFFERENT
    // chromosomes (multi-chrom paralog families, e.g. RABL2A/RABL2B) must not cluster together, and the
    // emitted event needs its chrom for the downstream microhomology check.
    let mut mosaics: Vec<(usize, usize, &str, u64, (u64, u64))> = calls
        .iter()
        .enumerate()
        .filter(|(_, c)| c.is_mosaic())
        .filter_map(|(i, c)| match (c.copy_a, c.copy_b, c.breakpoint_ref) {
            (Some(a), Some(b), Some(br)) => {
                Some((a, b, chroms.get(i).copied().unwrap_or(""), (br.0 + br.1) / 2, br))
            }
            _ => None,
        })
        .collect();
    // Stable order by (a, b, chrom, midpoint).
    mosaics.sort_by(|x, y| (x.0, x.1, x.2, x.3).cmp(&(y.0, y.1, y.2, y.3)));

    let mut events: Vec<ConversionEvent> = Vec::new();
    let mut i = 0usize;
    while i < mosaics.len() {
        let (a, b, chrom, _, _) = mosaics[i];
        // Greedily grow a cluster of same oriented pair ON THE SAME CHROM within breakpoint_tol.
        let mut j = i;
        let mut mids: Vec<u64> = Vec::new();
        let mut br_lo = u64::MAX;
        let mut br_hi = 0u64;
        while j < mosaics.len() {
            let (aj, bj, chromj, midj, brj) = mosaics[j];
            if aj != a || bj != b || chromj != chrom {
                break;
            }
            if let Some(&last) = mids.last() {
                if midj.saturating_sub(last) > p.breakpoint_tol {
                    break;
                }
            }
            mids.push(midj);
            br_lo = br_lo.min(brj.0);
            br_hi = br_hi.max(brj.1);
            j += 1;
        }
        let n = mids.len();
        let dispersion = mids.last().copied().unwrap_or(0).saturating_sub(mids[0]);
        let confirmed = n >= p.family_min_supporting_reads && dispersion <= p.max_breakpoint_dispersion;
        events.push(ConversionEvent {
            copy_a: a,
            copy_b: b,
            chrom: chrom.to_string(),
            breakpoint_ref: (br_lo, br_hi),
            n_supporting_reads: n,
            breakpoint_dispersion: dispersion,
            confirmed,
        });
        i = j;
    }
    events
}

/// Unified gene-conversion-vs-artifact verdict for a mosaic family event. `aggregate_family`'s
/// `confirmed` flag captures only ONE leg (recurrence across molecules); but recurrence alone is
/// insufficient — a sequence-driven template-switch hotspot (microhomology at the same point)
/// produces RT-switch chimeras that ALSO recur, so they pass the recurrence gate. The discriminator
/// therefore needs two ORTHOGONAL legs in addition to recurrence:
///   * **microhomology** at the breakpoint = the RT/template-switch signature (a direct repeat
///     flanking the switch point — `genome::is_rt_switch` applied to the breakpoint bracket);
///   * **DNA support** = heritability: a real (historical) gene conversion is in the genome and so
///     recurs in matched DNA reads; an RT/template switch is an RNA-library artifact, absent from DNA.
/// Both legs are passed in as `Option<bool>` so "no evidence available" (`None`) is distinct from
/// "negative evidence" (`Some(false)`). The DNA leg can act as a **veto** (`Some(false)` → `Ambiguous`)
/// when a RELIABLE absence source exists; this lets the two cheap legs (recurrence + microhomology)
/// ship without the DNA catalog wired, while a reliable DNA source strengthens or vetoes the call.
///
/// MEASURED (bench/mosaic_discriminator/dna_support.py, T2T DNA PSV catalog): the catalog signal is
/// real (42% of multi-copy families show a heritable-conversion DNA mosaic) but SPARSE and ref0-centric
/// (only ~2.9% of the genome is in a ref0 interval; localized mosaics return "absent" almost everywhere
/// even in families that HAVE one). So catalog "absent" is UNRELIABLE negative evidence — a catalog-
/// backed DNA closure must return `Some(true)` / `None` only (positive corroboration), NEVER
/// `Some(false)`, or it would wrongly downgrade real conversions. The production paths pass `None`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Classification {
    /// recurrent + present in DNA + no template-switch signature → real biological gene conversion.
    GeneConversion,
    /// direct-repeat/microhomology at the breakpoint and NOT DNA-confirmed → RT/template-switch artifact.
    RtSwitchArtifact,
    /// sporadic (did not recur) with no template-switch signature → one-off chimera, lean artifact.
    ChimeraSuspect,
    /// conflicting or insufficient evidence (e.g. microhomology AND DNA support; or recurrent but DNA unknown).
    Ambiguous,
}

/// Classify one family event from the three orthogonal legs (recurrence via `ev.confirmed`,
/// `microhomology`, `dna_supported`). Pure: the caller supplies the two genome/DNA-derived signals.
/// Evaluation order matters — the microhomology-artifact rule fires BEFORE the gene-conversion rule,
/// so a recurrent-but-RT-signature event is correctly called an artifact rather than a conversion.
pub fn classify_event(
    ev: &ConversionEvent,
    microhomology: Option<bool>,
    dna_supported: Option<bool>,
) -> Classification {
    let mh = microhomology == Some(true);
    let dna_present = dna_supported == Some(true);
    let dna_absent = dna_supported == Some(false);
    if mh && !dna_present {
        // template-switch signature, not rescued by positive DNA support → artifact (even if it recurs).
        Classification::RtSwitchArtifact
    } else if ev.confirmed && !mh && !dna_absent {
        // recurrent + no template signature + DNA not contradicting (present or unchecked) → conversion.
        Classification::GeneConversion
    } else if !ev.confirmed && !mh {
        // sporadic, no signature → one-off chimera.
        Classification::ChimeraSuspect
    } else {
        // microhomology∧DNA-present conflict, or recurrent-but-DNA-absent (heritability contradicted).
        Classification::Ambiguous
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build per-site obs from a label string where each char is the copy index that the read
    // uniquely matches; '*' = Ambig (all copies match), '.' = Novel (no copy matches).
    fn obs_from(labels: &str, n_copies: usize) -> Vec<SiteObs> {
        labels
            .chars()
            .enumerate()
            .map(|(i, ch)| {
                let mut match_bits = vec![false; n_copies];
                match ch {
                    '*' => match_bits.iter_mut().for_each(|m| *m = true),
                    '.' => {}
                    c => {
                        let idx = c.to_digit(10).unwrap() as usize;
                        match_bits[idx] = true;
                    }
                }
                SiteObs { ref_pos: 1000 + i as u64 * 10, match_bits }
            })
            .collect()
    }

    fn p() -> MosaicParams { MosaicParams::default() }

    #[test]
    fn clean_switch_is_mosaic() {
        let c = detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p());
        assert_eq!(c.status, MosaicStatus::Mosaic);
        assert_eq!((c.copy_a, c.copy_b), (Some(0), Some(1)));
        assert_eq!((c.tract_a_sites, c.tract_b_sites), (3, 3));
        assert_eq!(c.margin, 3);
        // breakpoint bracket between the last 0-site and the first 1-site.
        assert_eq!(c.breakpoint_ref, Some((1000 + 2 * 10, 1000 + 3 * 10)));
    }

    #[test]
    fn reversed_switch_records_direction() {
        let c = detect_mosaic(&obs_from("111000", 2), 2, 0.005, &p());
        assert!(c.is_mosaic());
        assert_eq!((c.copy_a, c.copy_b), (Some(1), Some(0)));
    }

    #[test]
    fn isolated_error_flip_is_not_mosaic() {
        // D=8, one interior flip; no a≠b contiguous split improves → margin 0.
        let c = detect_mosaic(&obs_from("00000100", 2), 2, 0.005, &p());
        assert!(!c.is_mosaic());
        assert!(c.margin < p().min_improvement);
    }

    #[test]
    fn pure_copy_is_single_copy() {
        let c = detect_mosaic(&obs_from("000000", 2), 2, 0.005, &p());
        assert_eq!(c.status, MosaicStatus::SingleCopy);
        assert!(!c.is_mosaic());
    }

    #[test]
    fn all_ambiguous_abstains_nonidentifiable() {
        let c = detect_mosaic(&obs_from("******", 2), 2, 0.005, &p());
        assert_eq!(c.status, MosaicStatus::NonIdentifiable);
    }

    #[test]
    fn too_few_sites_abstains_lowpower() {
        let c = detect_mosaic(&obs_from("0011", 2), 2, 0.005, &p());
        assert_eq!(c.status, MosaicStatus::LowPower);
    }

    #[test]
    fn short_right_tract_fails_min_tract() {
        // D=7, right tract only 2 sites (< min_tract_sites=3) → not Mosaic.
        let c = detect_mosaic(&obs_from("0000011", 2), 2, 0.005, &p());
        assert!(!c.is_mosaic());
        assert!(c.tract_b_sites < p().min_tract_sites);
    }

    #[test]
    fn novel_wildcard_does_not_break_tracts() {
        // '.' (Novel) between clean 3+3 tracts is neutral.
        let c = detect_mosaic(&obs_from("000.111", 2), 2, 0.005, &p());
        assert!(c.is_mosaic());
        assert_eq!((c.tract_a_sites, c.tract_b_sites), (3, 3));
    }

    #[test]
    fn ambiguous_sites_are_neutral_in_tracts() {
        let c = detect_mosaic(&obs_from("00*00111*1", 2), 2, 0.005, &p());
        assert!(c.is_mosaic());
        assert_eq!((c.copy_a, c.copy_b), (Some(0), Some(1)));
    }

    #[test]
    fn three_copy_family_picks_discriminating_pair() {
        // copies 0 and 2 are the switching pair in a 3-copy family.
        let c = detect_mosaic(&obs_from("000222", 3), 3, 0.005, &p());
        assert!(c.is_mosaic());
        assert_eq!((c.copy_a, c.copy_b), (Some(0), Some(2)));
    }

    #[test]
    fn de_none_uses_wider_eps_is_fail_safe() {
        // A borderline call must clear a HARDER bar at the larger (fail-safe) eps.
        let lo = detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p());
        let hi = detect_mosaic(&obs_from("000111", 2), 2, 0.05, &p());
        assert!(lo.score > hi.score); // larger eps → smaller LR → smaller margin-over-threshold
    }

    #[test]
    fn family_confirms_reproducible_breakpoint() {
        // 3 independent reads with the SAME switch/breakpoint → Confirmed.
        let calls: Vec<MosaicCall> = (0..3)
            .map(|_| detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p()))
            .collect();
        let events = aggregate_family(&calls, &vec!["c1"; calls.len()], &p());
        assert_eq!(events.len(), 1);
        assert!(events[0].confirmed);
        assert_eq!(events[0].n_supporting_reads, 3);
    }

    #[test]
    fn family_event_carries_its_chrom() {
        let calls: Vec<MosaicCall> = (0..3)
            .map(|_| detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p()))
            .collect();
        let events = aggregate_family(&calls, &vec!["chrX"; calls.len()], &p());
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].chrom, "chrX");
    }

    #[test]
    fn same_breakpoint_on_different_chroms_does_not_cluster() {
        // 3 reads on chrom A + 3 reads on chrom B, all with the IDENTICAL switch/midpoint. They must
        // form TWO separate events (one per chrom), not one merged cluster — the multi-chrom paralog
        // (e.g. RABL2A/RABL2B) case. Chrom is part of the cluster key.
        let calls: Vec<MosaicCall> = (0..6)
            .map(|_| detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p()))
            .collect();
        let chroms = ["cA", "cA", "cA", "cB", "cB", "cB"];
        let events = aggregate_family(&calls, &chroms, &p());
        assert_eq!(events.len(), 2, "different chroms must not cluster together");
        let mut cs: Vec<&str> = events.iter().map(|e| e.chrom.as_str()).collect();
        cs.sort();
        assert_eq!(cs, vec!["cA", "cB"]);
        assert!(events.iter().all(|e| e.confirmed && e.n_supporting_reads == 3));
    }

    #[test]
    fn family_rejects_singleton_as_chimera_suspect() {
        let calls = vec![detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p())];
        let events = aggregate_family(&calls, &vec!["c1"; calls.len()], &p());
        assert_eq!(events.len(), 1);
        assert!(!events[0].confirmed); // 1 molecule < family_min_supporting_reads
    }

    // ----- classify_event (unified gene-conversion-vs-artifact discriminator) -----

    fn ev(confirmed: bool) -> ConversionEvent {
        ConversionEvent {
            copy_a: 0,
            copy_b: 1,
            chrom: "c1".to_string(),
            breakpoint_ref: (1000, 1010),
            n_supporting_reads: if confirmed { 5 } else { 1 },
            breakpoint_dispersion: 0,
            confirmed,
        }
    }

    #[test]
    fn classify_recurrent_dna_no_microhomology_is_gene_conversion() {
        assert_eq!(
            classify_event(&ev(true), Some(false), Some(true)),
            Classification::GeneConversion
        );
    }

    #[test]
    fn classify_microhomology_without_dna_is_rt_switch_even_if_recurrent() {
        // the load-bearing case: recurrence ALONE would have called this confirmed, but the
        // template-switch signature (microhomology) + no DNA support overrides it to artifact.
        assert_eq!(
            classify_event(&ev(true), Some(true), Some(false)),
            Classification::RtSwitchArtifact
        );
        assert_eq!(
            classify_event(&ev(true), Some(true), None),
            Classification::RtSwitchArtifact
        );
    }

    #[test]
    fn classify_sporadic_no_signature_is_chimera_suspect() {
        assert_eq!(
            classify_event(&ev(false), Some(false), None),
            Classification::ChimeraSuspect
        );
    }

    #[test]
    fn classify_microhomology_and_dna_conflict_is_ambiguous() {
        // direct repeat AND present in DNA: could be a real conversion at a repeat-prone site — abstain.
        assert_eq!(
            classify_event(&ev(true), Some(true), Some(true)),
            Classification::Ambiguous
        );
    }

    #[test]
    fn classify_recurrent_no_microhomology_dna_unchecked_is_gene_conversion() {
        // DNA is a VETO, not a requirement: unchecked (None) does not block the two cheap legs.
        assert_eq!(
            classify_event(&ev(true), Some(false), None),
            Classification::GeneConversion
        );
    }

    #[test]
    fn classify_recurrent_but_dna_absent_is_ambiguous() {
        // DNA was CHECKED and the breakpoint is ABSENT from the genome → contradicts heritability.
        assert_eq!(
            classify_event(&ev(true), Some(false), Some(false)),
            Classification::Ambiguous
        );
    }

    // ----- ground-truth confusion matrix: full real path (detect -> aggregate -> genome
    //       microhomology -> classify), BOTH directions, over a constructed genome -----

    /// Build N recurrent recombinant reads that switch copy 0 -> copy 1 at site index `bp`, with
    /// their breakpoint bracket placed at genome positions `(left, right)` (so the genome's sequence
    /// at those coords decides microhomology). Returns the aggregated (single) confirmed event.
    fn recurrent_event_at(bp_left: u64, bp_right: u64, n_reads: usize) -> ConversionEvent {
        // 6 decisive sites in ASCENDING genomic order: three copy-0 sites ending exactly at `bp_left`,
        // then three copy-1 sites starting exactly at `bp_right` = a clean "000111" 0->1 switch whose
        // breakpoint bracket is (bp_left, bp_right).
        let positions = [
            (bp_left - 20, 0usize),
            (bp_left - 10, 0),
            (bp_left, 0),
            (bp_right, 1),
            (bp_right + 10, 1),
            (bp_right + 20, 1),
        ];
        let mut calls = Vec::new();
        for _ in 0..n_reads {
            let obs: Vec<SiteObs> = positions
                .iter()
                .map(|&(ref_pos, copy)| {
                    let mut mb = vec![false; 2];
                    mb[copy] = true;
                    SiteObs { ref_pos, match_bits: mb }
                })
                .collect();
            calls.push(detect_mosaic(&obs, 2, 0.005, &p()));
        }
        let mut events = aggregate_family(&calls, &vec!["c1"; calls.len()], &p());
        assert_eq!(events.len(), 1, "one oriented breakpoint cluster");
        let ev = events.pop().unwrap();
        assert!(ev.confirmed, "recurrent across {} molecules -> confirmed", n_reads);
        ev
    }

    #[test]
    fn ground_truth_conversion_vs_rt_switch_confusion_matrix() {
        use crate::genome::GenomeIndex;
        // One genome, two breakpoint loci:
        //  * CONVERSION locus: breakpoint flanks DIFFER -> no microhomology -> GeneConversion.
        //  * RT-SWITCH locus: breakpoint sits at an exact direct repeat -> microhomology -> artifact.
        let mut seq = vec![b'A'; 400];
        // RT-switch direct repeat: 8 bp ending at 200 == 8 bp ending at 240.
        seq[192..200].copy_from_slice(b"CGTACGTA");
        seq[232..240].copy_from_slice(b"CGTACGTA");
        // Conversion locus: flanks ending at 300 vs 340 deliberately DIFFER.
        seq[292..300].copy_from_slice(b"CGTACGTA");
        seq[332..340].copy_from_slice(b"TTGGAACC");
        let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);

        let mh = |left: u64, right: u64| g.breakpoint_microhomology("c1", left, right, 6, 12);

        // RT-switch breakpoint bracket (200, 240): direct repeat present.
        let rt_ev = recurrent_event_at(200, 240, 5);
        assert!(mh(200, 240), "direct repeat at the RT-switch breakpoint");
        assert_eq!(
            classify_event(&rt_ev, Some(mh(200, 240)), None),
            Classification::RtSwitchArtifact,
            "recurrent + microhomology + DNA-unchecked -> artifact (recurrence alone would have mis-called it)"
        );

        // Conversion breakpoint bracket (300, 340): no direct repeat.
        let gc_ev = recurrent_event_at(300, 340, 5);
        assert!(!mh(300, 340), "no direct repeat at the conversion breakpoint");
        assert_eq!(
            classify_event(&gc_ev, Some(mh(300, 340)), None),
            Classification::GeneConversion,
            "recurrent + no microhomology + DNA-not-contradicting -> gene conversion"
        );

        // DNA leg as a veto: same clean conversion, but ABSENT from DNA -> downgraded to Ambiguous.
        assert_eq!(
            classify_event(&gc_ev, Some(mh(300, 340)), Some(false)),
            Classification::Ambiguous,
        );
    }
}
