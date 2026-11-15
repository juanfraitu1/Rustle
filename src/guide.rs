//! Guided mode (register 1622, docs/PREREG_guided_assembly_2026-11-10.md, rule G1d): a guide GTF applied to a finished
//! de novo assembly. Guide chains with a full-chain primary alignment are added; de novo multi-exon chains outside the guide
//! are dropped when they are a fragment or an intron-retaining version of a guide chain, need NEW_INTRON_READS when they put
//! a new intron beside guide introns, and NOVEL_READS otherwise. Single-exon models and chains in the guide pass unchanged.

use std::collections::{HashMap, HashSet};

pub const NEW_INTRON_READS: u32 = 10;
pub const NOVEL_READS: u32 = 3;

/// Introns, 0-based half-open, sorted.
pub type Chain = Vec<(i64, i64)>;

pub fn chain(exons: &[(i64, i64)]) -> Chain {
    exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

/// The guide's multi-exon chains on one (chrom, strand), indexed by intron.
#[derive(Default)]
pub struct GuideSet {
    pub chains: HashSet<Chain>,
    by_intron: HashMap<(i64, i64), Vec<Chain>>,
}

impl GuideSet {
    pub fn insert(&mut self, ch: Chain) {
        if !ch.is_empty() && self.chains.insert(ch.clone()) {
            for i in &ch {
                self.by_intron.entry(*i).or_default().push(ch.clone());
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Verdict {
    /// single-exon, or its chain is a guide chain
    Keep,
    /// contiguous run of a guide chain
    Fragment,
    /// a guide chain with one or more introns left out
    Retained,
    /// new intron beside guide introns, reads below NEW_INTRON_READS
    NewIntronLow,
    /// outside the guide, reads below NOVEL_READS
    NovelLow,
}

/// `ch` taken from `g` in order (introns may be skipped); `contiguous`: no skipping.
fn within(ch: &Chain, g: &Chain, contiguous: bool) -> bool {
    if contiguous {
        return g.windows(ch.len()).any(|w| w == ch.as_slice());
    }
    let mut it = g.iter();
    ch.iter().all(|i| it.any(|j| j == i))
}

/// The verdict for a de novo model with this chain and read count; `guide` = the guide on the model's (chrom, strand).
pub fn verdict(guide: Option<&GuideSet>, ch: &Chain, reads: u32) -> Verdict {
    let floor = |r: u32| if r >= NOVEL_READS { Verdict::Keep } else { Verdict::NovelLow };
    if ch.is_empty() {
        return Verdict::Keep;
    }
    let Some(gs) = guide else { return floor(reads) };
    if gs.chains.contains(ch) {
        return Verdict::Keep;
    }
    let n = ch.iter().filter(|i| gs.by_intron.contains_key(i)).count();
    if n == ch.len() {
        // every guide chain that could contain `ch` holds its first intron
        let cand = &gs.by_intron[&ch[0]];
        if cand.iter().any(|g| within(ch, g, true)) {
            return Verdict::Fragment;
        }
        if cand.iter().any(|g| within(ch, g, false)) {
            return Verdict::Retained;
        }
        return floor(reads);
    }
    if n > 0 {
        return if reads >= NEW_INTRON_READS { Verdict::Keep } else { Verdict::NewIntronLow };
    }
    floor(reads)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(chains: &[&[(i64, i64)]]) -> GuideSet {
        let mut g = GuideSet::default();
        chains.iter().for_each(|c| g.insert(c.to_vec()));
        g
    }
    fn guide() -> GuideSet {
        set(&[&[(10, 20), (30, 40), (50, 60)], &[(10, 20), (50, 60)]])
    }

    #[test]
    fn guide_chains_and_single_exon_models_pass_at_any_support() {
        let g = guide();
        assert_eq!(verdict(Some(&g), &vec![(10, 20), (50, 60)], 1), Verdict::Keep);
        assert_eq!(verdict(Some(&g), &vec![], 1), Verdict::Keep);
    }

    #[test]
    fn a_contiguous_run_of_a_guide_chain_is_a_fragment_whatever_its_reads() {
        assert_eq!(verdict(Some(&guide()), &vec![(30, 40), (50, 60)], 500), Verdict::Fragment);
    }

    #[test]
    fn a_guide_chain_minus_an_intron_is_retained_unless_that_is_itself_a_guide_chain() {
        let g = set(&[&[(10, 20), (30, 40), (50, 60), (70, 80)]]);
        assert_eq!(verdict(Some(&g), &vec![(10, 20), (50, 60)], 500), Verdict::Retained);
        assert_eq!(verdict(Some(&g), &vec![(50, 60), (10, 20)], 500), Verdict::Keep); // not in the guide's order: novel
    }

    #[test]
    fn a_new_intron_beside_guide_introns_needs_ten_reads() {
        let g = guide();
        assert_eq!(verdict(Some(&g), &vec![(10, 20), (33, 40)], 9), Verdict::NewIntronLow);
        assert_eq!(verdict(Some(&g), &vec![(10, 20), (33, 40)], 10), Verdict::Keep);
    }

    #[test]
    fn chains_with_no_guide_intron_need_three_reads() {
        let g = guide();
        assert_eq!(verdict(Some(&g), &vec![(100, 200)], 2), Verdict::NovelLow);
        assert_eq!(verdict(Some(&g), &vec![(100, 200)], 3), Verdict::Keep);
        assert_eq!(verdict(None, &vec![(100, 200)], 2), Verdict::NovelLow);
    }

    #[test]
    fn chain_is_the_gaps_between_sorted_exons() {
        assert_eq!(chain(&[(0, 10), (20, 30), (40, 50)]), vec![(10, 20), (30, 40)]);
        assert!(chain(&[(0, 10)]).is_empty());
    }
}
