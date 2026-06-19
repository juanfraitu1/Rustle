//! Family-aware RESCUE thin-locus scan (integration stage 4b) — the BAM-neighbourhood orchestration of
//! `bench/family_rescue.py`.
//!
//! Recover an under-ASSEMBLED copy that falls below the general-assembly `>= 3`-read gate but is
//! sequence-homologous to a CONFIRMED de-novo family (borrow strength / partial pooling). The decision core
//! (`family_rescue::rescue_thin_locus`: canonical-k-mer pre-filter + both-orientation POA) is already
//! ported; this builds the thin candidate loci from reads in the family neighbourhood:
//!   1. group reads by exact intron chain, collapse OVERLAPPING chains into loci (best-supported chain wins);
//!   2. drop loci overlapping an already-assembled family member span;
//!   3. build each thin locus's spliced sequence (canonical-junction gate, reverse-complement for `-`);
//!   4. POA-confirm against the family members; dedup by locus, keep the best core coverage.

use std::collections::BTreeMap;

use super::denovo_assemble::{build_spliced_seq, PrimaryRead};
use super::family_rescue::{rescue_thin_locus, FamilyMember, RescueOutcome, RescueParams};
use crate::genome::GenomeIndex;

/// Minimum read support for a thin candidate locus (`family_rescue.py::MIN_SUPPORT`).
pub const RESCUE_MIN_SUPPORT: u32 = 1;
/// Minimum spliced length of a thin locus to attempt rescue.
pub const RESCUE_MIN_LEN: usize = 200;

/// A thin candidate locus: a (collapsed) intron chain with read support, below the general-assembly gate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ThinLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub support: u32,
    pub introns: Vec<(u64, u64)>,
}

/// An already-assembled family member's genomic span (used to EXCLUDE overlapping thin loci).
#[derive(Clone, Debug)]
pub struct MemberSpan {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

/// A rescued copy: the thin locus, the family/member that confirmed it (+ core coverage + orientation), and
/// the spliced sequence (so the rescued copy is assignable).
#[derive(Clone, Debug)]
pub struct RescuedCopy {
    pub locus: ThinLocus,
    pub outcome: RescueOutcome,
    pub seq: Vec<u8>,
}

/// Group reads by exact intron chain, then collapse OVERLAPPING chains (by span) into loci, keeping the
/// best-supported chain (by support, then span) as each locus's representative. Multi-exon only. Mirrors
/// `family_rescue.py`'s per-interval locus collapse. Deterministic (sorted by start, ties by chain).
pub fn thin_loci(reads: &[PrimaryRead], min_support: u32) -> Vec<ThinLocus> {
    // group reads by (chrom, exact intron chain) -> (support, min_start, max_end); multi-exon only.
    let mut chains: BTreeMap<(&str, Vec<(u64, u64)>), (u32, u64, u64)> = BTreeMap::new();
    for r in reads {
        if r.introns.is_empty() {
            continue;
        }
        let e = chains.entry((r.chrom.as_str(), r.introns.clone())).or_insert((0, u64::MAX, 0));
        e.0 += 1;
        e.1 = e.1.min(r.ref_start);
        e.2 = e.2.max(r.ref_end);
    }
    // per chrom: collapse overlapping chains (single-linkage by span); the rep is the best (support, span).
    let mut by_chrom: BTreeMap<&str, Vec<(Vec<(u64, u64)>, u32, u64, u64)>> = BTreeMap::new();
    for ((chrom, chain), (sup, s, e)) in chains {
        by_chrom.entry(chrom).or_default().push((chain, sup, s, e));
    }
    struct Loc {
        s: u64,
        e: u64,
        sup: u32,
        s2: u64,
        e2: u64,
        intr: Vec<(u64, u64)>,
    }
    let mut out = Vec::new();
    for (chrom, mut group) in by_chrom {
        group.sort_by_key(|&(_, _, s, _)| s); // stable: ties keep (chrom,chain) order
        let mut loci: Vec<Loc> = Vec::new();
        for (intr, sup, s, e) in group {
            if sup < min_support {
                continue;
            }
            let mut merged = false;
            for l in loci.iter_mut() {
                if s <= l.e && e >= l.s {
                    l.s = l.s.min(s);
                    l.e = l.e.max(e);
                    if (sup, e - s) > (l.sup, l.e2 - l.s2) {
                        l.sup = sup;
                        l.s2 = s;
                        l.e2 = e;
                        l.intr = intr.clone();
                    }
                    merged = true;
                    break;
                }
            }
            if !merged {
                loci.push(Loc { s, e, sup, s2: s, e2: e, intr });
            }
        }
        for l in loci {
            out.push(ThinLocus {
                chrom: chrom.to_string(),
                start: l.s2,
                end: l.e2,
                support: l.sup,
                introns: l.intr,
            });
        }
    }
    out
}

/// Rescue under-assembled copies: for each thin locus NOT overlapping a member span, build its spliced
/// sequence and POA-confirm it against the family `members` (`rescue_thin_locus`). Dedup by locus, keeping
/// the best `core_recip`. Mirrors `family_rescue.py`'s scan + dedup.
pub fn rescue_thin_loci(
    loci: &[ThinLocus],
    members: &[FamilyMember],
    member_spans: &[MemberSpan],
    genome: &GenomeIndex,
    p: &RescueParams,
) -> Vec<RescuedCopy> {
    let mut by_key: BTreeMap<(String, u64, u64), RescuedCopy> = BTreeMap::new();
    for locus in loci {
        // exclude a thin locus that overlaps an already-assembled family member span.
        if member_spans
            .iter()
            .any(|m| m.chrom == locus.chrom && locus.start < m.end && locus.end > m.start)
        {
            continue;
        }
        let (seq, _strand) = match build_spliced_seq(genome, &locus.chrom, locus.start, locus.end, &locus.introns) {
            Some(v) => v,
            None => continue,
        };
        if seq.len() < RESCUE_MIN_LEN || seq.len() > p.len_cap {
            continue;
        }
        if let Some(outcome) = rescue_thin_locus(&seq, members, p) {
            let key = (locus.chrom.clone(), locus.start, locus.end);
            let better = by_key
                .get(&key)
                .map_or(true, |prev| outcome.core_recip > prev.outcome.core_recip);
            if better {
                by_key.insert(key, RescuedCopy { locus: locus.clone(), outcome, seq });
            }
        }
    }
    by_key.into_values().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut rng = SplitMix64(seed);
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..n).map(|_| B[(rng.next_u64() % 4) as usize]).collect()
    }
    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.iter().flat_map(|p| p.iter().copied()).collect()
    }
    fn read(chrom: &str, s: u64, e: u64, introns: &[(u64, u64)]) -> PrimaryRead {
        PrimaryRead { chrom: chrom.into(), ref_start: s, ref_end: e, introns: introns.to_vec() }
    }

    /// Genome with one thin gene: exon1 [0,200), canonical intron [200,220) GT..AG, exon2 [220,420);
    /// spliced = flank(50) + `core`(300) + flank(50). Returns the genome.
    fn thin_gene_genome(core: &[u8]) -> GenomeIndex {
        let mut g = vec![b'A'; 500];
        g[0..50].copy_from_slice(&rand_seq(50, 0x71));
        g[50..200].copy_from_slice(&core[0..150]);
        g[200] = b'G';
        g[201] = b'T';
        g[218] = b'A';
        g[219] = b'G';
        g[220..370].copy_from_slice(&core[150..300]);
        g[370..420].copy_from_slice(&rand_seq(50, 0x72));
        GenomeIndex::from_seqs(&[("c1", &g)])
    }

    // ---- thin_loci ----

    #[test]
    fn thin_loci_groups_and_keeps_multi_exon() {
        let reads = [
            read("c1", 0, 400, &[(100, 200)]),
            read("c1", 5, 410, &[(100, 200)]),
            read("c1", 0, 400, &[]), // single-exon -> ignored
        ];
        let loci = thin_loci(&reads, 1);
        assert_eq!(loci.len(), 1);
        assert_eq!(loci[0].support, 2);
        assert_eq!(loci[0].introns, vec![(100, 200)]);
    }

    #[test]
    fn thin_loci_collapses_overlapping_chains_best_supported() {
        // chain A [0,400) support 2, chain B [100,500) support 1 overlap -> one locus, rep = chain A.
        let reads = [
            read("c1", 0, 400, &[(100, 200)]),
            read("c1", 0, 400, &[(100, 200)]),
            read("c1", 100, 500, &[(150, 250)]),
        ];
        let loci = thin_loci(&reads, 1);
        assert_eq!(loci.len(), 1, "overlapping chains collapse to one locus");
        assert_eq!(loci[0].support, 2, "best-supported chain wins");
        assert_eq!(loci[0].introns, vec![(100, 200)]);
        assert_eq!((loci[0].start, loci[0].end), (0, 400));
    }

    #[test]
    fn thin_loci_separates_nonoverlapping() {
        let reads = [
            read("c1", 0, 300, &[(100, 200)]),
            read("c1", 1000, 1300, &[(1100, 1200)]),
        ];
        assert_eq!(thin_loci(&reads, 1).len(), 2);
    }

    // ---- rescue_thin_loci ----

    #[test]
    fn rescue_recovers_thin_locus_homologous_to_family() {
        let core = rand_seq(300, 0xC0FE_F00D);
        let genome = thin_gene_genome(&core);
        // a family member sharing the same core (distinct flanks)
        let mseq = cat(&[&rand_seq(50, 0x81), &core, &rand_seq(50, 0x82)]);
        let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
        // a SINGLE read (support 1, below the >=3 gate) forming the thin locus
        let reads = [read("c1", 0, 420, &[(200, 220)])];
        let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
        assert_eq!(loci.len(), 1);
        let rescued = rescue_thin_loci(&loci, &members, &[], &genome, &RescueParams::default());
        assert_eq!(rescued.len(), 1, "the thin copy is rescued into the family");
        assert_eq!(rescued[0].outcome.family_id, "FAM1");
        assert!(rescued[0].outcome.core_recip >= 0.13);
        assert_eq!(rescued[0].seq.len(), 400);
    }

    #[test]
    fn rescue_excludes_loci_overlapping_a_member_span() {
        let core = rand_seq(300, 0xC0FE_F00D);
        let genome = thin_gene_genome(&core);
        let mseq = cat(&[&rand_seq(50, 0x81), &core, &rand_seq(50, 0x82)]);
        let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
        let reads = [read("c1", 0, 420, &[(200, 220)])];
        let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
        // a member already assembled across the locus span -> the thin locus is NOT a new copy.
        let spans = [MemberSpan { chrom: "c1".into(), start: 100, end: 300 }];
        let rescued = rescue_thin_loci(&loci, &members, &spans, &genome, &RescueParams::default());
        assert!(rescued.is_empty(), "locus overlapping an assembled member is excluded");
    }

    #[test]
    fn rescue_rejects_non_homologous_thin_locus() {
        let genome = thin_gene_genome(&rand_seq(300, 0xC0FE_F00D));
        // a family member with a DIFFERENT, unrelated core -> no rescue.
        let mseq = cat(&[&rand_seq(50, 0x81), &rand_seq(300, 0xDEAD_BEEF), &rand_seq(50, 0x82)]);
        let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
        let reads = [read("c1", 0, 420, &[(200, 220)])];
        let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
        let rescued = rescue_thin_loci(&loci, &members, &[], &genome, &RescueParams::default());
        assert!(rescued.is_empty(), "a non-homologous thin locus is not rescued");
    }
}
