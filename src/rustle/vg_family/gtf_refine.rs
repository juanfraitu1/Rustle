//! Opt-in refinement of `copy_assign --gtf`'s de novo isoform set (`--gtf-refine`).
//!
//! **STATUS:** OPT-IN
//!
//! Reached only from `src/bin/copy_assign.rs`'s `if args.gtf` block, behind `--gtf-refine` (default empty =
//! byte-identical). Four annotation-free components derived from the human chr20 error-pattern dissection;
//! rules and thresholds are frozen in `docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md`.
//! Not validated until the pre-registered held-out chr17 run.

use std::collections::HashMap;
#[allow(unused_imports)]
use std::collections::HashSet;

use crate::genome::GenomeIndex;
use crate::vg_family::denovo_assemble::{PrimaryRead, Skeleton};
#[allow(unused_imports)]
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
}
