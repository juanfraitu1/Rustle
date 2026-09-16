//! Opt-in refinement of `copy_assign --gtf`'s de novo isoform set (`--gtf-refine`).
//!
//! **STATUS:** OPT-IN
//!
//! Reached only from `src/bin/copy_assign.rs`'s `if args.gtf` block, behind `--gtf-refine` (default empty =
//! byte-identical). Four annotation-free components derived from the human chr20 error-pattern dissection;
//! rules and thresholds are frozen in `docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md`.
//! Not validated until the pre-registered held-out chr17 run.

use std::collections::{HashMap, HashSet};

use crate::genome::GenomeIndex;
use crate::vg_family::denovo_assemble::{PrimaryRead, Skeleton};
use crate::vg_family::family_detect::DenovoTranscript;

/// Strand of an intron chain when EVERY junction is canonical on one consistent strand, else `None` —
/// the strict assemble gate's rule (`build_spliced_seq`, strict branch).
pub fn strict_chain_strand(genome: &GenomeIndex, chrom: &str, introns: &[(u64, u64)]) -> Option<char> {
    let mut strand: Option<char> = None;
    for &(d, a) in introns {
        let js = if genome.is_canonical_junction(chrom, d, a, '+') {
            '+'
        } else if genome.is_canonical_junction(chrom, d, a, '-') {
            '-'
        } else {
            return None;
        };
        if strand.is_some_and(|s| s != js) {
            return None;
        }
        strand = Some(js);
    }
    strand
}

/// True iff `read` is a 3'-anchored, intron-compatible fragment of the chain `full` whose strand is
/// `full_strand`: a strictly shorter contiguous sub-chain, not reaching into `full`'s neighbouring introns,
/// ending at `full`'s 3' end ('+': last intron; '-': first intron).
fn is_three_prime_fragment(read: &PrimaryRead, full: &[(u64, u64)], full_strand: Option<char>) -> bool {
    let r = &read.introns;
    let m = r.len();
    let Some(st) = full_strand else { return false };
    if m == 0 || full.len() <= m {
        return false;
    }
    for i in 0..=(full.len() - m) {
        if full[i..i + m] != r[..] {
            continue;
        }
        if i > 0 && read.ref_start < full[i - 1].1 {
            continue;
        }
        if i + m < full.len() && read.ref_end > full[i + m].0 {
            continue;
        }
        if (st == '+' && i + m != full.len()) || (st == '-' && i != 0) {
            continue;
        }
        return true;
    }
    false
}

/// `fragsupport` (spec Part 2). `all_chains` must be `pass1_skeletons(reads, 1)`: every spliced chain with at
/// least one exact read, carrying pass-1's own boundary and strand fields. Returns the spliced, non-footprint
/// chains whose `exact + fragments >= min_support`, with `n_reads` set to that total. A read counts as a
/// fragment only when it is a fragment of EXACTLY ONE candidate (assign-or-abstain, no 1/k splitting);
/// candidates are looked up through the read's first intron. Unspliced reads never count.
pub fn fragment_supported_spliced(
    all_chains: &[Skeleton],
    reads: &[PrimaryRead],
    min_support: u32,
    chain_strand: impl Fn(&str, &[(u64, u64)]) -> Option<char>,
) -> Vec<Skeleton> {
    let cands: Vec<&Skeleton> = all_chains.iter().filter(|s| !s.introns.is_empty() && !s.footprint).collect();
    let mut by_intron: HashMap<(&str, (u64, u64)), Vec<usize>> = HashMap::new();
    for (ci, c) in cands.iter().enumerate() {
        for &it in &c.introns {
            by_intron.entry((c.chrom.as_str(), it)).or_default().push(ci);
        }
    }
    let strands: Vec<Option<char>> = cands.iter().map(|c| chain_strand(&c.chrom, &c.introns)).collect();
    let mut frag = vec![0u32; cands.len()];
    for r in reads {
        let Some(&first) = r.introns.first() else { continue };
        let Some(ids) = by_intron.get(&(r.chrom.as_str(), first)) else { continue };
        let hits: Vec<usize> =
            ids.iter().copied().filter(|&ci| is_three_prime_fragment(r, &cands[ci].introns, strands[ci])).collect();
        if hits.len() == 1 {
            frag[hits[0]] += 1;
        }
    }
    cands
        .iter()
        .enumerate()
        .filter_map(|(ci, c)| {
            let total = c.n_reads + frag[ci];
            (total >= min_support).then(|| Skeleton { n_reads: total, ..(*c).clone() })
        })
        .collect()
}

fn exons_of(start: u64, end: u64, introns: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut ex = Vec::with_capacity(introns.len() + 1);
    let mut prev = start;
    for &(d, a) in introns {
        ex.push((prev, d));
        prev = a;
    }
    ex.push((prev, end));
    ex
}

/// `subset` (spec Part 2): indices of multi-exon models whose intron chain is a contiguous proper sub-chain of
/// another model on the same chrom and strand, with both terminal exons inside the container's matching exons
/// (up to 5 bp of overhang allowed only into a container INTRON, never past the container's own ends).
/// Computed against the whole input at once; no support condition. (The proper-sub-chain requirement equals
/// the chr20 simulation whenever chains are unique, which the gate guarantees.)
pub fn subset_removals(models: &[DenovoTranscript]) -> HashSet<usize> {
    let mut idx: HashMap<(&str, char, (u64, u64)), Vec<(usize, usize)>> = HashMap::new();
    for (mi, m) in models.iter().enumerate() {
        for (p, &it) in m.introns.iter().enumerate() {
            idx.entry((m.chrom.as_str(), m.strand, it)).or_default().push((mi, p));
        }
    }
    let mut rm = HashSet::new();
    for (ti, t) in models.iter().enumerate() {
        let k = t.introns.len();
        if k == 0 {
            continue;
        }
        let Some(hits) = idx.get(&(t.chrom.as_str(), t.strand, t.introns[0])) else { continue };
        for &(si, p) in hits {
            let s = &models[si];
            if si == ti || s.introns.len() <= k || p + k > s.introns.len() || s.introns[p..p + k] != t.introns[..] {
                continue;
            }
            let sx = exons_of(s.start, s.end, &s.introns);
            let lo = sx[p].0 as i64 - t.start as i64;
            let ro = t.end as i64 - sx[p + k].1 as i64;
            let left_ok = lo <= 0 || (p != 0 && lo <= 5);
            let right_ok = ro <= 0 || (p + k != s.introns.len() && ro <= 5);
            if left_ok && right_ok {
                rm.insert(ti);
                break;
            }
        }
    }
    rm
}

/// `mono` rule 1 (spec Part 2): single-exon models lying fully inside an exon of a not-removed multi-exon
/// model on the same chrom and strand.
pub fn mono_in_own_exon_removals(models: &[DenovoTranscript], removed: &HashSet<usize>) -> HashSet<usize> {
    let mut exons_by: HashMap<(&str, char), Vec<(u64, u64)>> = HashMap::new();
    for (i, m) in models.iter().enumerate() {
        if removed.contains(&i) || m.introns.is_empty() {
            continue;
        }
        exons_by.entry((m.chrom.as_str(), m.strand)).or_default().extend(exons_of(m.start, m.end, &m.introns));
    }
    let mut rm = HashSet::new();
    for (i, m) in models.iter().enumerate() {
        if removed.contains(&i) || !m.introns.is_empty() {
            continue;
        }
        if let Some(ex) = exons_by.get(&(m.chrom.as_str(), m.strand)) {
            if ex.iter().any(|&(s, e)| s <= m.start && e >= m.end) {
                rm.insert(i);
            }
        }
    }
    rm
}

/// `mono` rule 2 (spec Part 2): single-exon model M (not already removed) is removed when spliced reads
/// dominate its locus: `SI >= U || SE >= U`, where over the reads overlapping M, `U` = unspliced reads,
/// `SI` = spliced reads with an intron covering >= 50% of M, `SE` = spliced reads on M's strand (FLAG 0x10)
/// with an aligned exon block covering >= 50% of M.
pub fn mono_spliced_dominated_removals(
    models: &[DenovoTranscript],
    removed: &HashSet<usize>,
    reads: &[PrimaryRead],
) -> HashSet<usize> {
    let mut rm = HashSet::new();
    for (i, m) in models.iter().enumerate() {
        if removed.contains(&i) || !m.introns.is_empty() {
            continue;
        }
        let len = (m.end - m.start) as i64;
        let ov = |a: u64, b: u64| b.min(m.end) as i64 - a.max(m.start) as i64;
        let (mut u, mut si, mut se) = (0usize, 0usize, 0usize);
        for r in reads {
            if r.chrom != m.chrom || r.ref_start >= m.end || r.ref_end <= m.start {
                continue;
            }
            if r.introns.is_empty() {
                u += 1;
                continue;
            }
            let intronic = r.introns.iter().map(|&(d, a)| ov(d, a)).max().unwrap_or(i64::MIN);
            if 2 * intronic >= len {
                si += 1;
            }
            let read_strand = if r.reverse { '-' } else { '+' };
            if read_strand == m.strand {
                let exonic = exons_of(r.ref_start, r.ref_end, &r.introns).iter().map(|&(a, b)| ov(a, b)).max().unwrap_or(i64::MIN);
                if 2 * exonic >= len {
                    se += 1;
                }
            }
        }
        if si >= u || se >= u {
            rm.insert(i);
        }
    }
    rm
}

/// Apply `subset`, then `mono` (rule 1 then rule 2) when enabled; returns the kept models in input order.
pub fn apply_post_filters(
    models: Vec<DenovoTranscript>,
    reads: &[PrimaryRead],
    subset: bool,
    mono: bool,
) -> Vec<DenovoTranscript> {
    let mut removed: HashSet<usize> = if subset { subset_removals(&models) } else { HashSet::new() };
    if mono {
        let r1 = mono_in_own_exon_removals(&models, &removed);
        removed.extend(r1);
        let r2 = mono_spliced_dominated_removals(&models, &removed, reads);
        removed.extend(r2);
    }
    models.into_iter().enumerate().filter(|(i, _)| !removed.contains(i)).map(|(_, m)| m).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(s: u64, e: u64, introns: Vec<(u64, u64)>) -> PrimaryRead {
        PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: e, introns, reverse: false }
    }

    fn chain(s: u64, e: u64, introns: Vec<(u64, u64)>, n: u32) -> Skeleton {
        Skeleton {
            chrom: "c1".into(), start: s, end: e, n_reads: n, introns, tied_seeded: false,
            read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: n,
        }
    }

    // full chain: exons [100,200) [300,400) [500,600) [700,800); introns (200,300) (400,500) (600,700)
    fn full() -> Skeleton {
        chain(100, 800, vec![(200, 300), (400, 500), (600, 700)], 1)
    }

    fn tx(s: u64, e: u64, strand: char, introns: Vec<(u64, u64)>) -> DenovoTranscript {
        DenovoTranscript { chrom: "c1".into(), start: s, end: e, strand, introns, ..Default::default() }
    }

    // container S: exons [100,200) [300,400) [500,600) [700,800)
    fn container() -> DenovoTranscript {
        tx(100, 800, '+', vec![(200, 300), (400, 500), (600, 700)])
    }

    #[test]
    fn plus_strand_suffix_fragment_rescues_a_one_read_chain() {
        let reads = vec![read(100, 800, full().introns.clone()), read(520, 800, vec![(600, 700)])];
        let out = fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+'));
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_reads, 2, "1 exact + 1 fragment");
        assert_eq!((out[0].start, out[0].end), (100, 800), "boundaries from exact reads only");
    }

    #[test]
    fn plus_strand_prefix_is_not_three_prime_anchored() {
        let reads = vec![read(100, 800, full().introns.clone()), read(100, 350, vec![(200, 300)])];
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).is_empty());
    }

    #[test]
    fn minus_strand_requires_the_prefix() {
        let reads = vec![read(100, 800, full().introns.clone()), read(100, 350, vec![(200, 300)])];
        assert_eq!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('-')).len(), 1);
    }

    #[test]
    fn a_read_extending_into_the_previous_intron_is_not_compatible() {
        // the read's chain (600,700) is full's suffix, but it starts at 380 < 500 = the acceptor of full's
        // preceding intron (400,500): it reaches into that intron, so it is not intron-compatible
        let reads = vec![read(100, 800, full().introns.clone()), read(380, 800, vec![(600, 700)])];
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).is_empty());
    }

    #[test]
    fn an_ambiguous_fragment_abstains() {
        let other = chain(550, 900, vec![(600, 700), (750, 850)], 1); // (600,700) is NOT its last intron
        let alt = chain(450, 800, vec![(500, 550), (600, 700)], 1); // (600,700) IS its last intron
        let reads = vec![
            read(100, 800, full().introns.clone()),
            read(450, 800, alt.introns.clone()),
            read(560, 800, vec![(600, 700)]), // suffix of BOTH full and alt -> abstain
        ];
        let out = fragment_supported_spliced(&[full(), alt, other], &reads, 2, |_, _| Some('+'));
        assert!(out.is_empty(), "no chain may take an ambiguous fragment");
    }

    #[test]
    fn unspliced_reads_and_unstranded_chains_give_no_support() {
        let reads = vec![read(100, 800, full().introns.clone()), read(720, 800, vec![])];
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).is_empty());
        let reads = vec![read(100, 800, full().introns.clone()), read(520, 800, vec![(600, 700)])];
        assert_eq!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).len(), 1, "control");
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| None).is_empty());
    }

    #[test]
    fn chains_already_at_support_are_kept_and_single_exon_or_footprint_candidates_are_ignored() {
        let two = chain(100, 800, vec![(200, 300), (400, 500), (600, 700)], 2);
        let mono = chain(1000, 1200, vec![], 5);
        let mut fp = chain(2000, 2600, vec![(2100, 2200)], 5);
        fp.footprint = true;
        let out = fragment_supported_spliced(&[two, mono, fp], &[], 2, |_, _| Some('+'));
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_reads, 2);
    }

    #[test]
    fn subset_removes_a_contained_sub_chain_and_keeps_the_container() {
        let t = tx(320, 800, '+', vec![(400, 500), (600, 700)]);
        assert_eq!(subset_removals(&[container(), t]), HashSet::from([1]));
    }

    #[test]
    fn subset_allows_5bp_overhang_only_into_a_container_intron() {
        // p = 1 (not the container's first exon): 5 bp into intron (200,300) is allowed, 6 is not
        assert_eq!(subset_removals(&[container(), tx(295, 800, '+', vec![(400, 500), (600, 700)])]), HashSet::from([1]));
        assert!(subset_removals(&[container(), tx(294, 800, '+', vec![(400, 500), (600, 700)])]).is_empty());
        // p = 0 (the container's first exon): any overhang past the container start is not allowed
        assert_eq!(subset_removals(&[container(), tx(100, 400, '+', vec![(200, 300)])]), HashSet::from([1]), "control");
        assert!(subset_removals(&[container(), tx(99, 400, '+', vec![(200, 300)])]).is_empty());
    }

    #[test]
    fn subset_never_crosses_strand_and_ignores_single_exon() {
        assert!(subset_removals(&[container(), tx(320, 800, '-', vec![(400, 500), (600, 700)])]).is_empty());
        assert!(subset_removals(&[container(), tx(320, 380, '+', vec![])]).is_empty());
    }

    #[test]
    fn mono_inside_own_spliced_exon_is_strand_aware() {
        let models = vec![container(), tx(510, 590, '+', vec![]), tx(510, 590, '-', vec![]), tx(150, 350, '+', vec![])];
        assert_eq!(mono_in_own_exon_removals(&models, &HashSet::new()), HashSet::from([1]));
        // a removed container no longer shelters anything
        assert!(mono_in_own_exon_removals(&models, &HashSet::from([0])).is_empty());
    }

    #[test]
    fn mono_spliced_dominance_counts_intronic_and_same_strand_exonic_reads() {
        let m = tx(420, 480, '+', vec![]); // len 60, inside intron (400,500)
        let spliced = PrimaryRead { chrom: "c1".into(), ref_start: 100, ref_end: 800, introns: container().introns.clone(), reverse: false };
        let unspliced = PrimaryRead { chrom: "c1".into(), ref_start: 410, ref_end: 490, introns: vec![], reverse: false };
        // 1 spliced read spans it with an intron vs 1 unspliced read -> SI (1) >= U (1) -> removed
        assert_eq!(mono_spliced_dominated_removals(&[m.clone()], &HashSet::new(), &[spliced.clone(), unspliced.clone()]), HashSet::from([0]));
        // 2 unspliced reads outvote it -> kept
        assert!(mono_spliced_dominated_removals(&[m.clone()], &HashSet::new(), &[spliced.clone(), unspliced.clone(), unspliced.clone()]).is_empty());
        // exonic evidence only counts on the model's strand
        let m2 = tx(520, 580, '+', vec![]); // inside exon [500,600)
        let rev = PrimaryRead { reverse: true, ..spliced.clone() };
        let u2 = PrimaryRead { chrom: "c1".into(), ref_start: 510, ref_end: 590, introns: vec![], reverse: false };
        assert!(mono_spliced_dominated_removals(&[m2.clone()], &HashSet::new(), &[rev, u2.clone()]).is_empty());
        assert_eq!(mono_spliced_dominated_removals(&[m2], &HashSet::new(), &[spliced, u2]), HashSet::from([0]));
    }

    #[test]
    fn apply_post_filters_runs_subset_then_mono_and_keeps_order() {
        let models = vec![container(), tx(320, 800, '+', vec![(400, 500), (600, 700)]), tx(510, 590, '+', vec![]), tx(5000, 5100, '+', vec![])];
        let u = PrimaryRead { chrom: "c1".into(), ref_start: 5000, ref_end: 5100, introns: vec![], reverse: false };
        let kept = apply_post_filters(models, &[u.clone(), u], true, true);
        let spans: Vec<(u64, u64)> = kept.iter().map(|t| (t.start, t.end)).collect();
        assert_eq!(spans, vec![(100, 800), (5000, 5100)]);
    }

    #[test]
    fn mono_with_no_overlapping_reads_is_removed() {
        let m = tx(420, 480, '+', vec![]);
        // Empty reads: U=SI=SE=0, so SI >= U (0 >= 0) is true -> removed (chr20 sim pins the >= operator)
        assert_eq!(mono_spliced_dominated_removals(&[m.clone()], &HashSet::new(), &[]), HashSet::from([0]));
        // Reads on different chromosome or not overlapping: same result (no overlaps -> U=0 -> removed)
        let non_overlap = PrimaryRead { chrom: "c2".into(), ref_start: 500, ref_end: 600, introns: vec![], reverse: false };
        assert_eq!(mono_spliced_dominated_removals(&[m], &HashSet::new(), &[non_overlap]), HashSet::from([0]));
    }

    #[test]
    fn subset_never_removes_either_of_two_identical_chains() {
        let dup = tx(150, 750, '+', container().introns.clone());
        assert!(subset_removals(&[container(), dup]).is_empty());
    }
}
