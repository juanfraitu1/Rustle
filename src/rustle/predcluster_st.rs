//! StringTie-faithful prediction SELECTION (sub-project 1). Default OFF
//! (RUSTLE_PREDCLUSTER_ST=1). Runs ST's selection sub-stages in ST's order on the
//! candidate predictions Rustle's flow already produced. See
//! docs/superpowers/specs/2026-05-30-predcluster-selection-parity-design.md.
use crate::bpcov::Bpcov;
use crate::path_extract::Transcript;
use crate::types::RunConfig;

// ── StringTie constants (rlink.h / stringtie.cpp) ────────────────────────────
/// `readthr` (stringtie.cpp:144): minimum per-bundle-bp read coverage to accept a
/// prediction. Default 1.0 in `-L` mode.
const READTHR: f64 = 1.0;
/// `mintranscriptlen` (stringtie.cpp:139): minimum exonic length to print.
const MINTRANSCRIPTLEN: u64 = 200;
/// `DROP` (rlink.h:13): included-prediction coverage drop factor.
const DROP: f64 = 0.5;
/// `ERROR_PERC` (rlink.h:14): error-tolerance fraction.
const ERROR_PERC: f64 = 0.1;
/// `CHI_WIN` (rlink.h:17): coverage window size; also the high-cov isofrac threshold.
const CHI_WIN: f64 = 100.0;
/// `CHI_THR` (rlink.h:18): coverage threshold; mid-cov isofrac relaxation point.
const CHI_THR: f64 = 50.0;
/// `isofrac` (stringtie.cpp:163): per-maxint long-read coverage fraction.
const ISOFRAC: f64 = 0.01;

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

/// Faithful port of ST's `retainedintron(pred, n1, n2, lowintron)` (rlink.cpp:17117).
///
/// `n1` is the container/killer (higher cov), `n2` the candidate victim (contained).
/// `low_n1` is `n1`'s per-intron low-coverage mask (from `build_lowintron_flags`,
/// itself the faithful port of ST's lowintron construction at rlink.cpp:18131-18167).
///
/// Returns true (kill) when, for some intron `i-1` of `n1` that is low-covered,
/// an exon of `n2` bridges that intron region under ST's coverage conditions:
///   - last exon of n2 overlapping the intron start, with cov < frac*cov_n1 → kill (1)
///   - a middle exon of n2 inside the intron → kill (2)
///   - first/edge exon overlap with cov < frac*cov_n1 → kill (1)
///
/// `frac` is `ERROR_PERC` in long-read mode (mixedMode isofrac relaxation does not
/// apply for `-L`). This is the actual gate the task asked us to determine: the RI
/// kill is **lowintron-gated** (per-base coverage), not unconditional.
fn retainedintron_st(
    n1: &Transcript,
    n2: &Transcript,
    low_n1: &[bool],
    frac: f64,
) -> bool {
    if n1.exons.len() < 2 || n2.exons.is_empty() {
        return false;
    }
    let a = &n1.exons; // killer exons
    let b = &n2.exons; // victim exons
    let mut j = 0usize;
    for i in 1..a.len() {
        if j > b.len().saturating_sub(1) {
            return false;
        }
        if !low_n1.get(i - 1).copied().unwrap_or(false) {
            continue;
        }
        // ST: if(j==n2->exons.Count()-1 && cov_n2<frac*cov_n1 && n2->exons[j].start<=n1->exons[i-1].end)
        if j == b.len() - 1 && n2.coverage < frac * n1.coverage && b[j].0 <= a[i - 1].1 {
            return true;
        }
        // ST: while(j<count && n2->exons[j].end<n1->exons[i].start) j++;
        while j < b.len() && b[j].1 < a[i].0 {
            j += 1;
        }
        // ST: if(!j && cov_n2<frac*cov_n1) return 1;
        if j == 0 && n2.coverage < frac * n1.coverage {
            return true;
        }
        // ST: if(j<count && n2->exons[j].start<=n1->exons[i-1].end)
        if j < b.len() && b[j].0 <= a[i - 1].1 {
            if j > 0 && j < b.len() - 1 {
                return true; // middle exon bridges the intron
            } else if n2.coverage < frac * n1.coverage {
                return true;
            }
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

/// Emit a `pred_kill` parity event matching ST's predcluster payload shape for the
/// pairwise stage (rlink.cpp:18537/18597): `reason`, `cov`, `nexons`, `stage:"pairwise"`.
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

/// Emit an isofrac-stage `pred_kill` event matching ST's payload at rlink.cpp:18782:
/// `reason:"isofrac"`, `cov`, `usedcov`, `multicov`, `isofraclong`, `nexons`,
/// `stage:"isofrac"`.
fn emit_isofrac_kill(t: &Transcript, usedcov: f64, multicov: f64, isofraclong: f64) {
    if !crate::parity::decisions::is_enabled() {
        return;
    }
    let pk_s = t.exons.first().map(|(s, _)| *s + 1).unwrap_or(0);
    let pk_e = t.exons.last().map(|(_, e)| *e).unwrap_or(0);
    let payload = format!(
        r#""reason":"isofrac","cov":{:.4},"usedcov":{:.4},"multicov":{:.4},"isofraclong":{:.6},"nexons":{},"stage":"isofrac""#,
        t.coverage,
        usedcov,
        multicov,
        isofraclong,
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
/// Near-equal/exact chain collapse (Task 5) is added later.
pub fn select_predictions_st(
    candidates: Vec<Transcript>,
    bpcov: Option<&Bpcov>,
    config: &RunConfig,
) -> Vec<Transcript> {
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

    // ── lowintron flags (ST: rlink.cpp:18131-18167) ──────────────────────────
    // Per-prediction, per-intron low-coverage mask, reproduced from bpcov via the
    // faithful port in transcript_filter::build_lowintron_flags. This is what gates
    // ST's retained-intron kill (the RI kill is NOT unconditional — see the design
    // note on retainedintron_st). Without bpcov we get an empty mask → RI never
    // fires, which is the correct degenerate behaviour (no per-base signal).
    let longreads = config.long_reads;
    let lowintron: Vec<Vec<bool>> = match bpcov {
        Some(bp) => crate::transcript_filter::build_lowintron_flags(
            &survivors,
            bp,
            longreads,
            config.singlethr,
            config.pairwise_error_perc,
            config.pairwise_drop,
        ),
        None => survivors
            .iter()
            .map(|t| vec![false; t.exons.len().saturating_sub(1)])
            .collect(),
    };

    // ── (B) pairwise containment (lowintron-gated RI + included_drop) ────────
    // Sort by coverage DESC (ST predord: highest coverage first). Stable so equal
    // covs keep input order, matching ST's tie behaviour closely enough for parity.
    let mut order: Vec<usize> = (0..survivors.len()).collect();
    order.sort_by(|&a, &b| {
        survivors[b]
            .coverage
            .partial_cmp(&survivors[a].coverage)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

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
            if tj.chrom != ti.chrom {
                continue;
            }
            // ST gates the whole pairwise block on `overlaps.get(n1,n2)`, a
            // *significant*-overlap test (substantial exonic overlap, not mere
            // span touching). spans_overlap is a coarse proxy; included_pred and
            // retainedintron_st below impose the real geometric structure.
            if !spans_overlap(tj, ti) {
                continue;
            }
            // ST pairwise branch chain (same/cross strand, n2 = ti the lower-cov):
            //   else if(retainedintron(n1,n2,lowintron)) → retained_intron
            //   ... (cross-strand handled separately) ...
            //   else if(exons_n2<=exons_n1 && cov_n1>cov_n2*DROP
            //           && included_pred(n1,n2)) → included_drop
            // The branches are mutually exclusive (`else if`). RI is checked first
            // and does NOT require structural inclusion — only that one of j's
            // introns is low-covered and i bridges it (retainedintron_st). Strand
            // must match for both kills (ST's RI is inside the same-strand path via
            // the lowintron geometry; the cross-strand branch has separate rules we
            // do not port here, so we require equal strand).
            if tj.strand != ti.strand {
                continue;
            }
            if retainedintron_st(tj, ti, &lowintron[j], ERROR_PERC) {
                dropped[i] = true;
                emit_pred_kill(ti, "retained_intron");
            } else if ti.exons.len() <= tj.exons.len()
                && tj.coverage > ti.coverage * DROP
                && crate::transcript_filter::included_pred(
                    &survivors,
                    j,
                    i,
                    longreads,
                    bpcov,
                    config.singlethr,
                    config.pairwise_error_perc,
                    false,
                )
            {
                dropped[i] = true;
                emit_pred_kill(ti, "included_drop");
            }
        }
    }

    // ── (C) per-maxint isofrac / longunder (ST: rlink.cpp:18734-18800) ───────
    // Long-read branch: per maximal-coverage interval, seed usedcov/multicov from
    // the dominant (highest abs(tlen)*cov) flagged prediction, then kill any later
    // prediction whose cov falls below isofraclong * usedcov (or multicov). Guides
    // (t_eq) are never killed. Only runs in long-read mode, matching ST's
    // `if(longreads) { ... }` gate.
    if longreads {
        isofrac_st(&survivors, &mut dropped);
    }

    survivors
        .into_iter()
        .enumerate()
        .filter(|(idx, _)| !dropped[*idx])
        .map(|(_, t)| t)
        .collect()
}

/// Maximal-coverage interval over a set of predictions: the elementary segments
/// `[bp[k], bp[k+1])` cut at every distinct exon boundary, each carrying the list of
/// predictions whose exons overlap it. ST builds these incrementally via
/// `add_exon_to_maxint` (rlink.cpp:17370) producing maximal runs; splitting at every
/// boundary is a strict refinement with identical per-interval membership, so the
/// per-interval dominant/longunder decisions are unchanged (kills are idempotent).
struct MaxInt {
    /// indices into the prediction slice whose exons overlap this interval.
    members: Vec<usize>,
}

fn build_maxint(txs: &[Transcript], alive: &[bool]) -> Vec<MaxInt> {
    let mut bps: Vec<u64> = Vec::new();
    for (n, t) in txs.iter().enumerate() {
        if !alive[n] {
            continue;
        }
        for &(s, e) in &t.exons {
            if e > s {
                bps.push(s);
                bps.push(e);
            }
        }
    }
    if bps.len() < 2 {
        return Vec::new();
    }
    bps.sort_unstable();
    bps.dedup();
    let nseg = bps.len() - 1;
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); nseg];
    for (n, t) in txs.iter().enumerate() {
        if !alive[n] {
            continue;
        }
        for &(es, ee) in &t.exons {
            if ee <= es {
                continue;
            }
            let si = bps.partition_point(|&x| x < es);
            let ei = bps.partition_point(|&x| x < ee);
            for seg in members.iter_mut().take(ei.min(nseg)).skip(si) {
                // avoid duplicate membership when a pred has two exons in one seg
                if seg.last() != Some(&n) {
                    seg.push(n);
                }
            }
        }
    }
    members
        .into_iter()
        .filter(|m| !m.is_empty())
        .map(|members| MaxInt { members })
        .collect()
}

/// Per-maxint isofrac / longunder kill, faithful to rlink.cpp:18734-18800 (longreads
/// branch). Marks `dropped[n]=true` for predictions that fall under the dominant
/// prediction's coverage by more than `isofraclong`.
fn isofrac_st(survivors: &[Transcript], dropped: &mut [bool]) {
    // membership only needs the still-alive preds.
    let alive: Vec<bool> = dropped.iter().map(|&d| !d).collect();
    let maxints = build_maxint(survivors, &alive);
    for mi in &maxints {
        // ST exord: flagged preds touching this interval, sorted by abs(tlen)*cov desc.
        let mut exord: Vec<usize> = mi
            .members
            .iter()
            .copied()
            .filter(|&n| !dropped[n])
            .collect();
        if exord.is_empty() {
            continue;
        }
        let priority = |n: usize| -> f64 {
            (exonic_len(&survivors[n]) as f64) * survivors[n].coverage
        };
        exord.sort_by(|&a, &b| {
            priority(b)
                .partial_cmp(&priority(a))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        // ST usedcov/multicov are indexed by strand [0]=-, [1]=+. Seed from exord[0].
        let mut usedcov = [0.0f64; 2];
        let mut multicov = [0.0f64; 2];
        let dom = exord[0];
        let s0 = if survivors[dom].strand == '+' { 1 } else { 0 };
        usedcov[s0] = survivors[dom].coverage;
        if survivors[dom].exons.len() > 1 {
            multicov[s0] = usedcov[s0];
        }
        for &n in exord.iter().skip(1) {
            if dropped[n] {
                continue;
            }
            let t = &survivors[n];
            let s = if t.strand == '+' { 1 } else { 0 };
            let cov = t.coverage;
            let mut longunder = false;
            let mut isofraclong = ISOFRAC;
            if !is_guide(t) {
                if cov > CHI_WIN {
                    isofraclong = ISOFRAC * ERROR_PERC * DROP;
                } else if cov > CHI_THR {
                    isofraclong = ISOFRAC * DROP;
                }
                if t.exons.len() > 1 {
                    if (multicov[s] == 0.0
                        && cov < isofraclong * usedcov[s]
                        && cov < DROP / ERROR_PERC)
                        || cov < isofraclong * multicov[s]
                    {
                        longunder = true;
                    }
                } else if cov < isofraclong * usedcov[s] {
                    longunder = true;
                }
            }
            if longunder {
                dropped[n] = true;
                emit_isofrac_kill(t, usedcov[s], multicov[s], isofraclong);
            } else {
                usedcov[s] += cov;
                if t.exons.len() > 1 {
                    multicov[s] += cov;
                }
            }
        }
    }
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

    /// Long-read config (-L) for select_predictions_st tests.
    fn cfg() -> RunConfig {
        let mut c = RunConfig::default();
        c.long_reads = true;
        c
    }

    /// Pass-through wrapper: no bpcov (lowintron empty → RI never fires).
    fn select(cands: Vec<Transcript>) -> Vec<Transcript> {
        select_predictions_st(cands, None, &cfg())
    }

    #[test]
    fn survival_gate_drops_low_cov() {
        // 3-exon, exonic len 300 > MINTRANSCRIPTLEN, cov 0.5 < READTHR → dropped.
        let cands = vec![tx('+', &[(0, 100), (200, 300), (400, 500)], 0.5)];
        assert_eq!(select(cands).len(), 0);
    }

    #[test]
    fn survival_gate_drops_short() {
        // cov ok but exonic len 50 < MINTRANSCRIPTLEN → dropped.
        let cands = vec![tx('+', &[(0, 50)], 10.0)];
        assert_eq!(select(cands).len(), 0);
    }

    #[test]
    fn survival_gate_keeps_guide_low_cov() {
        let cands = vec![guide('+', &[(0, 50)], 0.0)];
        assert_eq!(select(cands).len(), 1);
    }

    #[test]
    fn retained_intron_gated_on_lowintron() {
        // ST's RI kill is lowintron-gated: with a low-coverage mask on j's first
        // intron (100,200), the contained i (which bridges it) is killed. cov_i 2.0
        // < ERROR_PERC*cov_j (2.0) is false → not killed via the cov branch, but the
        // middle-exon branch (j>0 && j<count-1) kills unconditionally when low.
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 20.0);
        let i = tx('+', &[(0, 300), (400, 500)], 0.5);
        // intron 0 of j = (100,200) is low; intron 1 = (300,400) not low.
        let low_j = vec![true, false];
        assert!(retainedintron_st(&j, &i, &low_j, ERROR_PERC));
        // Without a low intron, RI does not fire.
        let low_none = vec![false, false];
        assert!(!retainedintron_st(&j, &i, &low_none, ERROR_PERC));
    }

    #[test]
    fn included_drop_unconditional_containment() {
        // contained i: exact sub-chain of j (introns ⊆, no retained intron, no low
        // mask). cov_j (10) > cov_i (8) * DROP (0.5) = 4 → included_drop.
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 10.0);
        let i = tx('+', &[(0, 100), (200, 300)], 8.0);
        let kept = select(vec![j, i]);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].exons.len(), 3);
    }

    #[test]
    fn different_strand_not_contained() {
        let j = tx('+', &[(0, 100), (200, 300), (400, 500)], 10.0);
        let i = tx('-', &[(0, 100), (200, 300)], 8.0);
        assert_eq!(select(vec![j, i]).len(), 2);
    }

    #[test]
    fn isofrac_kills_underwater_alt_variant() {
        // Two multi-exon predictions that overlap in a maxint interval but have
        // NON-subset intron chains (so pairwise containment does NOT fire and the
        // kill must come from isofrac). They share exon [0,100) (first maxint) but
        // diverge afterwards: dom splices (100,200); alt splices (100,250).
        // Dom cov 4000 (>CHI_WIN → isofraclong=0.0005). In [0,100) dom dominant,
        // multicov=4000. alt cov 1.05 (≥READTHR). multicov branch:
        // 1.05 < 0.0005*4000 = 2.0 → longunder → alt killed.
        let dom = tx('+', &[(0, 100), (200, 600)], 4000.0);
        let alt = tx('+', &[(0, 100), (250, 600)], 1.05);
        // not subset (alt's intron (100,250) ∉ dom) → pairwise keeps both.
        let kept = select(vec![dom, alt]);
        assert_eq!(kept.len(), 1);
        assert!((kept[0].coverage - 4000.0).abs() < 1e-9);
    }

    #[test]
    fn isofrac_keeps_comparable_cov() {
        // Same overlap geometry, comparable coverage. Dom cov 40 (<CHI_THR →
        // isofraclong=ISOFRAC=0.01); multicov=40; alt cov 20 < 0.01*40 = 0.4? no →
        // both kept (and pairwise can't fire: non-subset chains).
        let dom = tx('+', &[(0, 100), (200, 600)], 40.0);
        let alt = tx('+', &[(0, 100), (250, 600)], 20.0);
        let kept = select(vec![dom, alt]);
        assert_eq!(kept.len(), 2);
    }

    #[test]
    fn guides_exempt_from_isofrac() {
        // A guide with tiny cov overlapping a high-cov dominant's maxint is never
        // killed by isofrac (ST: !pred[n]->t_eq gate). Non-subset chains so pairwise
        // also can't drop it.
        let dom = tx('+', &[(0, 100), (200, 600)], 4000.0);
        let g = guide('+', &[(0, 100), (250, 600)], 0.5);
        let kept = select(vec![dom, g]);
        assert!(kept.iter().any(|t| t.ref_transcript_id.is_some()));
    }
}
