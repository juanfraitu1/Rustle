//! StringTie-faithful prediction SELECTION (sub-project 1). Default OFF
//! (RUSTLE_PREDCLUSTER_ST=1). Runs ST's selection sub-stages in ST's order on the
//! candidate predictions Rustle's flow already produced. See
//! docs/superpowers/specs/2026-05-30-predcluster-selection-parity-design.md.
use crate::path_extract::Transcript;

// ── StringTie constants (rlink.h / stringtie.cpp) ────────────────────────────
/// `readthr` (stringtie.cpp:144): minimum per-bundle-bp read coverage to accept a
/// prediction. Default 1.0 in `-L` mode.
const READTHR: f64 = 1.0;
/// `mintranscriptlen` (stringtie.cpp:139): minimum exonic length to print.
const MINTRANSCRIPTLEN: u64 = 200;
/// `DROP` (rlink.h:13): included-prediction coverage drop factor.
const DROP: f64 = 0.5;

/// Half-open exonic length of a transcript (Σ exon lengths).
fn exonic_len(t: &Transcript) -> u64 {
    t.exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum()
}

/// Guide-ness mirrors `pred->t_eq` (ST) / `is_guide_pair` (Rustle): a transcript is
/// a guide when it carries a reference id or a `guide:` source tag. Guides are never
/// dropped by the survival gate or the pairwise containment filter.
fn is_guide(t: &Transcript) -> bool {
    t.ref_transcript_id.is_some()
        || t
            .source
            .as_deref()
            .map_or(false, |s| s.starts_with("guide:"))
}

/// Intron-coordinate set of a transcript: `(exons[k].1, exons[k+1].0)` for each
/// adjacent exon pair. ST compares intron *patterns* via node bitvecs; deriving the
/// intron coordinates from `exons` is equivalent for same-graph predictions and
/// avoids threading graph node ids into this pure function (per the design doc).
fn intron_chain(t: &Transcript) -> Vec<(u64, u64)> {
    t.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

/// `i`'s intron set ⊆ `j`'s intron set (ST: `(p_i.b & p_j.b) == p_i.b`).
/// Single-exon `i` (empty intron set) is trivially contained.
fn introns_subset(i_introns: &[(u64, u64)], j_introns: &[(u64, u64)]) -> bool {
    i_introns.iter().all(|x| j_introns.contains(x))
}

/// Faithful port of ST's `has_retained_intron(p1=j, p2=i)` (rlink.cpp:7191): returns
/// true when container `j` has an intron edge that contained `i` lacks but spans —
/// i.e. `i` retained an intron that `j` splices out. We detect this structurally: an
/// intron `(d, a)` in `j` for which `i` has no matching intron, but `i` has an exon
/// that overlaps the intron interval `(d, a)` (the retained-intron exon bridges it).
fn has_retained_intron(j: &Transcript, i: &Transcript) -> bool {
    let j_introns = intron_chain(j);
    let i_introns = intron_chain(i);
    for &(d, a) in &j_introns {
        if i_introns.contains(&(d, a)) {
            continue; // i also splices this intron — not retained
        }
        // Does any exon of i span across j's intron (d, a)? An exon (es, ee) spans
        // the intron when it covers strictly inside (d, a): ee > d && es < a.
        if i.exons.iter().any(|&(es, ee)| ee > d && es < a) {
            return true;
        }
    }
    false
}

/// Two transcripts overlap on the genome (same chrom assumed; checked by caller).
fn spans_overlap(a: &Transcript, b: &Transcript) -> bool {
    let (as_, ae) = match (a.exons.first(), a.exons.last()) {
        (Some(f), Some(l)) => (f.0, l.1),
        _ => return false,
    };
    let (bs, be) = match (b.exons.first(), b.exons.last()) {
        (Some(f), Some(l)) => (f.0, l.1),
        _ => return false,
    };
    as_ <= be && bs <= ae
}

/// Emit a `pred_kill` parity event matching ST's predcluster payload shape
/// (rlink.cpp:18537/18597): `reason`, `cov`, `nexons`, `stage:"pairwise"`.
fn emit_pred_kill(t: &Transcript, reason: &str) {
    if !crate::parity::decisions::is_enabled() {
        return;
    }
    let pk_s = t.exons.first().map(|(s, _)| *s + 1).unwrap_or(0);
    let pk_e = t.exons.last().map(|(_, e)| *e).unwrap_or(0);
    let payload = format!(
        r#""reason":"{}","cov":{:.4},"nexons":{},"stage":"pairwise""#,
        reason,
        t.coverage,
        t.exons.len()
    );
    crate::parity::decisions::emit("pred_kill", Some(&t.chrom), pk_s, pk_e, t.strand, &payload);
}

/// StringTie-faithful prediction selection.
///
/// Runs ST's selection sub-stages, in ST's order, on the candidate predictions
/// Rustle's flow already produced. This task implements the first two:
///
/// (A) cov/readthr/length survival gate (ST: rlink.cpp:449/465 + the `readthr`
///     check on `pred->cov`); guides are exempt.
/// (B) pairwise containment in coverage-descending order, mirroring ST's
///     `print_predcluster` pairwise loop (rlink.cpp:18494-18613):
///       - container `j` (higher cov), contained `i` (lower cov), same strand,
///         overlapping, `i`'s intron set ⊆ `j`'s intron set, `i` not a guide;
///       - if `has_retained_intron(j, i) && cov_i < cov_j` → drop i (retained_intron);
///       - else if `cov_j > cov_i * DROP` (the included_pred/DROP branch) → drop i
///         (included_drop). This covers the unconditional-containment case (the
///         container is virtually always > 0.5× the contained cov).
///
/// Isofrac (Task 4) and near-equal/exact chain collapse (Task 5) are added later.
pub fn select_predictions_st(candidates: Vec<Transcript>) -> Vec<Transcript> {
    // ── (A) cov/readthr/length survival gate ─────────────────────────────────
    let survivors: Vec<Transcript> = candidates
        .into_iter()
        .filter(|t| {
            if is_guide(t) {
                return true; // guides are never dropped by the survival gate
            }
            if t.coverage < READTHR {
                return false;
            }
            if exonic_len(t) < MINTRANSCRIPTLEN {
                return false;
            }
            true
        })
        .collect();

    // ── (B) pairwise containment (RI unconditional + included_drop) ──────────
    // Sort by coverage DESC (ST predord: highest coverage first). Stable so equal
    // covs keep input order, matching ST's tie behaviour closely enough for parity.
    let mut order: Vec<usize> = (0..survivors.len()).collect();
    order.sort_by(|&a, &b| {
        survivors[b]
            .coverage
            .partial_cmp(&survivors[a].coverage)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Precompute intron chains once.
    let chains: Vec<Vec<(u64, u64)>> = survivors.iter().map(intron_chain).collect();

    let mut dropped = vec![false; survivors.len()];

    // Outer: container candidate j (higher cov). Inner: contained candidate i
    // (lower cov, appears later in coverage-descending order).
    for oi in 0..order.len() {
        let j = order[oi];
        if dropped[j] {
            continue;
        }
        for &i in order.iter().skip(oi + 1) {
            if dropped[i] {
                continue;
            }
            // i must be the contained (lower-priority) prediction and not a guide.
            if is_guide(&survivors[i]) {
                continue;
            }
            let (tj, ti) = (&survivors[j], &survivors[i]);
            if tj.strand != ti.strand || tj.chrom != ti.chrom {
                continue;
            }
            if !spans_overlap(tj, ti) {
                continue;
            }
            // i's exon count must be <= j's (ST: pred[n2]->exons.Count()<=pred[n1]->...).
            if ti.exons.len() > tj.exons.len() {
                continue;
            }
            // i's intron set ⊆ j's intron set.
            if !introns_subset(&chains[i], &chains[j]) {
                continue;
            }
            // ST pairwise branch order: retained_intron first, then included_drop.
            if has_retained_intron(tj, ti) {
                if ti.coverage < tj.coverage {
                    dropped[i] = true;
                    emit_pred_kill(ti, "retained_intron");
                }
                // RI present but cov_i >= cov_j: ST keeps i (no kill).
            } else if tj.coverage > ti.coverage * DROP {
                // Contained, no retained intron → ST's included_pred/DROP kill.
                dropped[i] = true;
                emit_pred_kill(ti, "included_drop");
            }
        }
    }

    survivors
        .into_iter()
        .enumerate()
        .filter(|(idx, _)| !dropped[*idx])
        .map(|(_, t)| t)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(strand: char, exons: &[(u64, u64)], cov: f64) -> Transcript {
        Transcript {
            chrom: "chr1".into(),
            strand,
            exons: exons.to_vec(),
            coverage: cov,
            ..Default::default()
        }
    }

    fn guide(strand: char, exons: &[(u64, u64)], cov: f64) -> Transcript {
        let mut t = tx(strand, exons, cov);
        t.ref_transcript_id = Some("REF1".into());
        t
    }

    #[test]
    fn survival_gate_drops_low_cov() {
        // 3-exon, exonic len 300 > MINTRANSCRIPTLEN, cov 0.5 < READTHR → dropped.
        let cands = vec![tx('+', &[(0, 100), (200, 300), (400, 500)], 0.5)];
        assert_eq!(select_predictions_st(cands).len(), 0);
    }

    #[test]
    fn survival_gate_drops_short() {
        // cov ok but exonic len 50 < MINTRANSCRIPTLEN → dropped.
        let cands = vec![tx('+', &[(0, 50)], 10.0)];
        assert_eq!(select_predictions_st(cands).len(), 0);
    }

    #[test]
    fn survival_gate_keeps_guide_low_cov() {
        let cands = vec![guide('+', &[(0, 50)], 0.0)];
        assert_eq!(select_predictions_st(cands).len(), 1);
    }

    #[test]
    fn retained_intron_dropped_when_lower_cov() {
        // container j: 3 exons, introns (100,200) and (300,400).
        // contained i: 2 exons; intron set {(300,400)} ⊆ j, and i's exon (0,300)
        // bridges j's (100,200) intron → retained intron. cov_i < cov_j → dropped.
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 20.0);
        let i = tx('+', &[(0, 300), (400, 500)], 2.0);
        let kept = select_predictions_st(vec![j, i]);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].coverage, 20.0);
    }

    #[test]
    fn included_drop_unconditional_containment() {
        // contained i: exact sub-chain of j (introns ⊆, no retained intron).
        // cov_j (10) > cov_i (8) * DROP (0.5) = 4 → included_drop.
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 10.0);
        let i = tx('+', &[(0, 100), (200, 300)], 8.0);
        let kept = select_predictions_st(vec![j, i]);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].exons.len(), 3);
    }

    #[test]
    fn different_strand_not_contained() {
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 10.0);
        let i = tx('-', &[(0, 100), (200, 300)], 8.0);
        assert_eq!(select_predictions_st(vec![j, i]).len(), 2);
    }

    #[test]
    fn ri_kept_when_higher_cov() {
        // Retained-intron variant i has HIGHER cov than the spliced j. Sorted desc,
        // i is the container-position; j (5) is tested as contained against i — but
        // j is NOT a subset of i (j has intron (100,200) that i lacks), so no kill.
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 5.0);
        let i = tx('+', &[(0, 300), (400, 500)], 20.0);
        assert_eq!(select_predictions_st(vec![j, i]).len(), 2);
    }
}
