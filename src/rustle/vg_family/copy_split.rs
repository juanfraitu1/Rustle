//! Joint read-coherence + PSV "copy split" decomposition.
//!
//! Two orthogonal axes group a family's reads into emitted (copy, isoform) units:
//!   1. STRUCTURAL axis (read-coherence): the exact ordered intron chain a read traverses.
//!      Different chains are different transcripts regardless of allele content.
//!   2. COPY axis (PSV haplotype): within one chain, reads that span the family's PSV columns
//!      carry the allele base distinguishing the paralog copy they came from. When >= 2 distinct
//!      allele groups differ at >= K columns, the chain SPLITS into multiple identifiable copies.
//!
//! The split is deliberately conservative against the indistinguishability wall: a lone sporadic
//! single-column mismatch never spawns a copy, copies separated by < K columns merge (identifiable
//! = false), and fully-ambiguous (all-None) reads never force a copy of their own.

use std::collections::BTreeMap;

/// One read's observation for joint copy+isoform grouping.
#[derive(Clone, Debug)]
pub struct ReadObs {
    /// read-coherence STRUCTURAL key: ordered intron chain (donor, acceptor).
    pub intron_chain: Vec<(u64, u64)>,
    /// COPY axis: allele base at each family PSV column (index = column; None if the read
    /// does not span that column). length == n_psv_columns.
    pub psv_alleles: Vec<Option<u8>>,
}

/// One emitted (copy, isoform): an intron chain + the PSV haplotype that distinguishes the copy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CopyIsoform {
    pub intron_chain: Vec<(u64, u64)>,
    pub allele_vector: Vec<Option<u8>>, // consensus PSV haplotype of this copy (None = unobserved/merged)
    pub read_count: usize,
    pub identifiable: bool,             // true if produced by a PSV split; false if a merged (non-identifiable) chain group
}

/// Joint read-coherence + PSV decomposition.
/// 1. group reads by EXACT intron_chain (read-coherence).
/// 2. within a chain-group: among reads that span the PSVs, form candidate copies = distinct
///    allele vectors with >= min_reads_per_copy supporting reads. If >= 2 candidate copies that
///    PAIRWISE differ at >= min_psv_k columns exist, SPLIT: emit one CopyIsoform per such copy
///    (identifiable=true); reads not matching a candidate copy (sub-threshold vectors, sporadic
///    single-column disagreements, or fully-ambiguous all-None reads) are apportioned to the
///    consistent copy if unique, else left shared (counted toward the group but not forcing a copy).
/// 3. otherwise emit ONE merged CopyIsoform for the chain-group (identifiable=false), read_count=all.
/// Deterministic: sort output by (intron_chain, allele_vector). Use BTreeMap, no HashMap-into-output.
pub fn split_readchain_by_psv(
    reads: &[ReadObs],
    n_psv_columns: usize,
    min_psv_k: usize,            // K: identifiability threshold (distinct PSV columns)
    min_reads_per_copy: usize,
) -> Vec<CopyIsoform> {
    // 1. Group reads by EXACT intron chain (read-coherence / structural axis).
    //    BTreeMap keeps the structural axis deterministically ordered.
    let mut by_chain: BTreeMap<Vec<(u64, u64)>, Vec<&ReadObs>> = BTreeMap::new();
    for r in reads {
        by_chain.entry(r.intron_chain.clone()).or_default().push(r);
    }

    let mut out: Vec<CopyIsoform> = Vec::new();

    for (chain, group) in by_chain {
        // 2. Within a chain-group, tally distinct fully-observed-at-their-columns allele
        //    vectors among reads that span at least one PSV column. A candidate copy is a
        //    distinct allele vector with >= min_reads_per_copy supporting reads. BTreeMap
        //    keeps copy enumeration deterministic.
        let mut vector_counts: BTreeMap<Vec<Option<u8>>, usize> = BTreeMap::new();
        for r in &group {
            // Skip fully-ambiguous reads: they span no column and can never define a copy.
            if r.psv_alleles.iter().all(|a| a.is_none()) {
                continue;
            }
            *vector_counts.entry(r.psv_alleles.clone()).or_insert(0) += 1;
        }

        // Candidate copies: distinct vectors with enough support.
        let candidates: Vec<Vec<Option<u8>>> = vector_counts
            .iter()
            .filter(|(_, &c)| c >= min_reads_per_copy)
            .map(|(v, _)| v.clone())
            .collect();

        // A valid identifiable split needs the PSV axis to actually carry >= K columns,
        // and >= 2 candidate copies that PAIRWISE differ at >= K columns.
        let split_copies: Vec<Vec<Option<u8>>> = if n_psv_columns >= min_psv_k {
            select_identifiable_copies(&candidates, min_psv_k)
        } else {
            Vec::new()
        };

        if split_copies.len() >= 2 {
            // 2a. SPLIT. Apportion every read in the group: a read joins a copy iff that copy
            //     is the UNIQUE copy consistent with the read's observed (non-None) columns.
            //     Reads consistent with >1 copy (e.g. all-None, or sub-threshold ambiguous
            //     vectors) are left shared and do not force or inflate any copy's count.
            let mut counts = vec![0usize; split_copies.len()];
            for r in &group {
                let mut consistent: Option<usize> = None;
                let mut unique = true;
                for (i, copy) in split_copies.iter().enumerate() {
                    if read_consistent_with(&r.psv_alleles, copy) {
                        if consistent.is_some() {
                            unique = false;
                            break;
                        }
                        consistent = Some(i);
                    }
                }
                if unique {
                    if let Some(i) = consistent {
                        counts[i] += 1;
                    }
                }
            }

            for (copy, count) in split_copies.into_iter().zip(counts) {
                out.push(CopyIsoform {
                    intron_chain: chain.clone(),
                    allele_vector: copy,
                    read_count: count,
                    identifiable: true,
                });
            }
        } else {
            // 3. MERGED chain-group: not identifiable. read_count = all reads in the group.
            //    allele_vector = consensus haplotype (per-column majority; None where unobserved
            //    or no majority emerges), purely informational since the copy isn't split out.
            let allele_vector = consensus_haplotype(&group, n_psv_columns);
            out.push(CopyIsoform {
                intron_chain: chain.clone(),
                allele_vector,
                read_count: group.len(),
                identifiable: false,
            });
        }
    }

    // Deterministic ordering: (intron_chain, allele_vector). by_chain already ordered chains;
    // within a chain copies were enumerated from a BTreeMap (allele-vector ordered). A final
    // stable sort makes the contract explicit regardless of construction order.
    out.sort_by(|a, b| {
        a.intron_chain
            .cmp(&b.intron_chain)
            .then_with(|| a.allele_vector.cmp(&b.allele_vector))
    });
    out
}

// ===========================================================================
// ReadObs BRIDGE: build a ReadObs from a read's spliced alignment + the family's
// PSV genomic positions. The load-bearing primitive is `allele_at`, which reads a
// read's base at a REFERENCE position by walking the CIGAR (the error-prone part).
// ===========================================================================

/// One read's spliced alignment (minimal model): ref_start (0-based), CIGAR ops as (op,len)
/// with op in {'M','I','D','N','S'} (match/ins/del/intron/softclip; '=','X' treated as M),
/// and the read sequence (no hard-clipped bases).
#[derive(Clone, Debug)]
pub struct AlignedRead {
    pub ref_start: u64,
    pub cigar: Vec<(char, u64)>,
    pub seq: Vec<u8>,
}

/// Read base aligned to reference position ref_pos (0-based), or None if ref_pos is not a
/// matched position in this read (inside an intron N, a deletion D, or outside the read span).
pub fn allele_at(read: &AlignedRead, ref_pos: u64) -> Option<u8> {
    // Walk the CIGAR, tracking the current reference coordinate and read (seq) coordinate.
    // Only M/=/X ops consume BOTH and align a seq base to a ref position; we return that
    // base when the ref coordinate hits ref_pos. N (intron) and D (deletion) advance the
    // ref coordinate but consume no read base -> any ref_pos inside them is unmatched (None).
    // I (insertion) and S (softclip) advance the read coordinate but consume no ref.
    let mut ref_cur = read.ref_start;
    let mut seq_cur: u64 = 0;
    for &(op, len) in &read.cigar {
        match op {
            'M' | '=' | 'X' => {
                // ref_pos in [ref_cur, ref_cur+len) maps to seq[seq_cur + (ref_pos-ref_cur)].
                if ref_pos >= ref_cur && ref_pos < ref_cur + len {
                    let off = ref_pos - ref_cur;
                    return read.seq.get((seq_cur + off) as usize).copied();
                }
                ref_cur += len;
                seq_cur += len;
            }
            'N' | 'D' => {
                // Consumes reference only; no read base aligns here.
                if ref_pos >= ref_cur && ref_pos < ref_cur + len {
                    return None;
                }
                ref_cur += len;
            }
            'I' | 'S' => {
                // Consumes read only; reference coordinate unchanged.
                seq_cur += len;
            }
            _ => {
                // Unknown op (e.g. 'H' hard-clip, 'P' pad): consume nothing here. Hard-clipped
                // bases are not present in seq, and pads touch neither coordinate.
            }
        }
    }
    None
}

/// Intron chain = the (ref_end_of_block, ref_start_of_next_block) gaps produced by N ops.
pub fn intron_chain_of(read: &AlignedRead) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut ref_cur = read.ref_start;
    for &(op, len) in &read.cigar {
        match op {
            'M' | '=' | 'X' | 'D' => {
                // Reference-consuming, alignment-contiguous ops (a D does not break the chain).
                ref_cur += len;
            }
            'N' => {
                // Intron: gap from current ref position to current+len.
                out.push((ref_cur, ref_cur + len));
                ref_cur += len;
            }
            'I' | 'S' => {
                // Read-only ops do not move the reference coordinate.
            }
            _ => {}
        }
    }
    out
}

/// Build a ReadObs: intron chain from N ops; psv_alleles[i] = allele_at(read, psv_positions[i]).
pub fn build_read_obs(read: &AlignedRead, psv_positions: &[u64]) -> ReadObs {
    ReadObs {
        intron_chain: intron_chain_of(read),
        psv_alleles: psv_positions.iter().map(|&p| allele_at(read, p)).collect(),
    }
}

/// True if a read's observed (non-None) PSV columns all agree with `copy`'s alleles.
/// A None in the read means "did not span" -> imposes no constraint at that column.
fn read_consistent_with(read: &[Option<u8>], copy: &[Option<u8>]) -> bool {
    read.iter().zip(copy.iter()).all(|(r, c)| match (r, c) {
        (Some(rb), Some(cb)) => rb == cb,
        (Some(_), None) => false, // read observes a base where the copy haplotype is unobserved
        (None, _) => true,        // read did not span this column: no constraint
    })
}

/// Greedily select a maximal set of candidate copies that PAIRWISE differ at >= K observed
/// columns. Candidates are pre-sorted (BTreeMap order); a new candidate is admitted only if it
/// is >= K-distinct from every already-admitted copy, guaranteeing the emitted copies are
/// mutually identifiable. Sub-K-distinct near-duplicates fall to the indistinguishability wall.
fn select_identifiable_copies(
    candidates: &[Vec<Option<u8>>],
    min_psv_k: usize,
) -> Vec<Vec<Option<u8>>> {
    let mut chosen: Vec<Vec<Option<u8>>> = Vec::new();
    for cand in candidates {
        if chosen
            .iter()
            .all(|c| distinct_columns(cand, c) >= min_psv_k)
        {
            chosen.push(cand.clone());
        }
    }
    chosen
}

/// Count columns where both vectors observe a base and the bases differ.
fn distinct_columns(a: &[Option<u8>], b: &[Option<u8>]) -> usize {
    a.iter()
        .zip(b.iter())
        .filter(|(x, y)| match (x, y) {
            (Some(xb), Some(yb)) => xb != yb,
            _ => false,
        })
        .count()
}

/// Per-column majority allele over a merged chain-group; None where no read observes the column
/// or no strict majority exists. Deterministic (BTreeMap tally, lowest base wins ties).
fn consensus_haplotype(group: &[&ReadObs], n_psv_columns: usize) -> Vec<Option<u8>> {
    let mut out = Vec::with_capacity(n_psv_columns);
    for col in 0..n_psv_columns {
        let mut tally: BTreeMap<u8, usize> = BTreeMap::new();
        for r in group {
            if let Some(Some(base)) = r.psv_alleles.get(col) {
                *tally.entry(*base).or_insert(0) += 1;
            }
        }
        // pick the base with the highest count; ties broken by lowest base (BTreeMap order).
        let best = tally
            .iter()
            .max_by(|x, y| x.1.cmp(y.1).then_with(|| y.0.cmp(x.0)))
            .map(|(b, _)| *b);
        out.push(best);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(chain: &[(u64, u64)], alleles: &[Option<u8>]) -> ReadObs {
        ReadObs {
            intron_chain: chain.to_vec(),
            psv_alleles: alleles.to_vec(),
        }
    }

    // Allele shorthands.
    const A: Option<u8> = Some(b'A');
    const C: Option<u8> = Some(b'C');
    const G: Option<u8> = Some(b'G');
    const T: Option<u8> = Some(b'T');
    const N: Option<u8> = None;

    /// Two reference intron chains to play the structural axis.
    fn chain_x() -> Vec<(u64, u64)> {
        vec![(100, 200), (300, 400)]
    }
    fn chain_y() -> Vec<(u64, u64)> {
        vec![(100, 200), (300, 450)]
    }

    #[test]
    fn split_two_copies_same_chain() {
        let cx = chain_x();
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(obs(&cx, &[A, C, T]));
        }
        for _ in 0..3 {
            reads.push(obs(&cx, &[G, T, A]));
        }
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 2, "expected exactly 2 copies");
        assert!(out.iter().all(|c| c.identifiable), "both must be identifiable");
        let counts: Vec<usize> = out.iter().map(|c| c.read_count).collect();
        assert_eq!(counts, vec![3, 3], "read counts 3 and 3");
        assert_ne!(out[0].allele_vector, out[1].allele_vector, "distinct allele vectors");
        assert!(out.iter().all(|c| c.intron_chain == cx));
    }

    #[test]
    fn no_oversplit_identical_alleles() {
        let cx = chain_x();
        let reads: Vec<ReadObs> = (0..6).map(|_| obs(&cx, &[A, C, T])).collect();
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 1, "single copy, no over-split");
        assert_eq!(out[0].read_count, 6);
        assert!(!out[0].identifiable, "single copy is not identifiable");
    }

    #[test]
    fn different_chains_are_distinct() {
        let cx = chain_x();
        let cy = chain_y();
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(obs(&cx, &[A, C, T]));
        }
        for _ in 0..3 {
            reads.push(obs(&cy, &[A, C, T]));
        }
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 2, "two distinct transcripts keyed by chain");
        let chains: Vec<&Vec<(u64, u64)>> = out.iter().map(|c| &c.intron_chain).collect();
        assert!(chains.contains(&&cx));
        assert!(chains.contains(&&cy));
        // identical alleles must NOT merge across chains.
        assert_ne!(out[0].intron_chain, out[1].intron_chain);
    }

    #[test]
    fn sporadic_error_does_not_split() {
        let cx = chain_x();
        let mut reads: Vec<ReadObs> = (0..5).map(|_| obs(&cx, &[A, C, T])).collect();
        reads.push(obs(&cx, &[A, C, A])); // lone variant differs at last col only
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 1, "lone variant neither >=K-distinct nor >=min_reads");
        assert_eq!(out[0].read_count, 6);
    }

    #[test]
    fn below_k_columns_merges() {
        let cx = chain_x();
        let mut reads: Vec<ReadObs> = (0..3).map(|_| obs(&cx, &[A, C, T])).collect();
        reads.extend((0..3).map(|_| obs(&cx, &[A, C, A]))); // differ at exactly 1 column
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 1, "differ at < K columns -> merge");
        assert!(!out[0].identifiable, "the indistinguishability wall");
    }

    #[test]
    fn ambiguous_read_not_forced() {
        let cx = chain_x();
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(obs(&cx, &[A, C, T]));
        }
        for _ in 0..3 {
            reads.push(obs(&cx, &[G, T, A]));
        }
        for _ in 0..2 {
            reads.push(obs(&cx, &[N, N, N])); // fully ambiguous
        }
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 2, "all-None reads create no third copy");
        assert!(out.iter().all(|c| c.identifiable));
    }

    // ---- Adversarial cases added by review ----

    /// Indistinguishability wall: zero PSV columns => no allele information at all =>
    /// exactly one merged, non-identifiable chain group regardless of read count.
    #[test]
    fn adv_no_psv_columns_merges_unidentifiable() {
        let cx = chain_x();
        // Reads carry empty allele vectors (n_psv_columns == 0).
        let reads: Vec<ReadObs> = (0..8).map(|_| obs(&cx, &[])).collect();
        let out = split_readchain_by_psv(&reads, 0, 2, 2);
        assert_eq!(out.len(), 1, "no PSV info => single merged copy");
        assert!(!out[0].identifiable, "n_psv=0 can never be identifiable");
        assert_eq!(out[0].read_count, 8, "all reads counted in the merged group");
        assert!(out[0].allele_vector.is_empty(), "no columns => empty haplotype");
    }

    /// A minority allele group with < min_reads_per_copy support must NOT become its own
    /// copy. Here [G,T,A] has only 1 supporting read (< min_reads=2) while [A,C,T] has 3.
    /// Result must be a single merged group, and the lone read folded into its count.
    #[test]
    fn adv_minority_below_min_reads_not_a_copy() {
        let cx = chain_x();
        let mut reads: Vec<ReadObs> = (0..3).map(|_| obs(&cx, &[A, C, T])).collect();
        reads.push(obs(&cx, &[G, T, A])); // only 1 read for a >=K-distinct vector
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 1, "lone strong-but-undersupported variant is not a copy");
        assert!(!out[0].identifiable, "only one candidate copy => not identifiable");
        assert_eq!(out[0].read_count, 4, "minority read folded into the merged count");
    }

    /// Apportionment: a partially-observed read consistent with exactly one copy is
    /// assigned to it; a fully-ambiguous read is left shared (no copy inflation).
    #[test]
    fn adv_partial_read_apportioned_uniquely() {
        let cx = chain_x();
        let mut reads = Vec::new();
        for _ in 0..2 {
            reads.push(obs(&cx, &[A, C, T]));
        }
        for _ in 0..2 {
            reads.push(obs(&cx, &[G, T, A]));
        }
        // [A,N,N] disagrees with [G,T,A] at col0 -> consistent ONLY with [A,C,T].
        reads.push(obs(&cx, &[A, N, N]));
        // [N,N,N] consistent with both -> shared, counts toward neither copy.
        reads.push(obs(&cx, &[N, N, N]));
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 2, "two identifiable copies");
        assert!(out.iter().all(|c| c.identifiable));
        let total: usize = out.iter().map(|c| c.read_count).sum();
        // 2 + 2 + 1 apportioned = 5; the all-None read is shared and excluded.
        assert_eq!(total, 5, "ambiguous read must not inflate any copy");
        // The copy whose col0==A must have picked up the [A,N,N] read (3), the other 2.
        let mut counts: Vec<usize> = out.iter().map(|c| c.read_count).collect();
        counts.sort();
        assert_eq!(counts, vec![2, 3], "partial read joined the unique consistent copy");
    }

    /// Determinism / order-independence: shuffling the input read order yields an
    /// identical output vector (same order, same counts).
    #[test]
    fn adv_deterministic_under_reorder() {
        let cx = chain_x();
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(obs(&cx, &[A, A, A]));
        }
        for _ in 0..3 {
            reads.push(obs(&cx, &[C, C, C]));
        }
        for _ in 0..3 {
            reads.push(obs(&cx, &[G, G, G]));
        }
        let out1 = split_readchain_by_psv(&reads, 3, 2, 2);
        reads.reverse();
        let out2 = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out1, out2, "output independent of input read order");
    }

    #[test]
    fn three_copies() {
        let cx = chain_x();
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(obs(&cx, &[A, A, A]));
        }
        for _ in 0..3 {
            reads.push(obs(&cx, &[C, C, C]));
        }
        for _ in 0..3 {
            reads.push(obs(&cx, &[G, G, G]));
        }
        let out = split_readchain_by_psv(&reads, 3, 2, 2);
        assert_eq!(out.len(), 3, "three pairwise >=K-distinct copies");
        assert!(out.iter().all(|c| c.identifiable));
    }

    // ---- ReadObs BRIDGE tests (CIGAR-walking allele_at is the load-bearing primitive) ----

    fn aligned(ref_start: u64, cigar: &[(char, u64)], seq: &[u8]) -> AlignedRead {
        AlignedRead {
            ref_start,
            cigar: cigar.to_vec(),
            seq: seq.to_vec(),
        }
    }

    /// "10M": straight match. Interior ref position returns the right base; before the
    /// start and past the end return None.
    #[test]
    fn bridge_allele_at_simple_match() {
        // ref 100..110 <-> seq[0..10]. seq = ACGTACGTAC
        let read = aligned(100, &[('M', 10)], b"ACGTACGTAC");
        assert_eq!(allele_at(&read, 100), Some(b'A'), "first matched ref pos -> seq[0]");
        assert_eq!(allele_at(&read, 103), Some(b'T'), "interior pos ref 103 -> seq[3]");
        assert_eq!(allele_at(&read, 109), Some(b'C'), "last matched pos -> seq[9]");
        assert_eq!(allele_at(&read, 99), None, "before ref_start -> None");
        assert_eq!(allele_at(&read, 110), None, "past the end -> None");
        assert_eq!(allele_at(&read, 200), None, "far past the end -> None");
    }

    /// "5M100N5M": a spliced read. Positions inside the N intron return None; positions in
    /// the second block return the correct base (ref advanced past the intron, read did NOT).
    #[test]
    fn bridge_allele_at_spliced_intron() {
        // block1: ref 100..105 <-> seq[0..5] = "AAAAA"
        // intron: ref 105..205 (N100), consumes no read
        // block2: ref 205..210 <-> seq[5..10] = "CCCGT"
        let read = aligned(100, &[('M', 5), ('N', 100), ('M', 5)], b"AAAAACCCGT");
        // exon1
        assert_eq!(allele_at(&read, 100), Some(b'A'));
        assert_eq!(allele_at(&read, 104), Some(b'A'), "last base of block1 -> seq[4]");
        // inside intron
        assert_eq!(allele_at(&read, 105), None, "intron start -> None");
        assert_eq!(allele_at(&read, 150), None, "deep in intron -> None");
        assert_eq!(allele_at(&read, 204), None, "last intron pos -> None");
        // exon2: ref 205 -> seq[5], NOT seq[105]; read coord did not advance over the intron.
        assert_eq!(allele_at(&read, 205), Some(b'C'), "block2 start -> seq[5]");
        assert_eq!(allele_at(&read, 208), Some(b'G'), "block2 interior ref 208 -> seq[8]");
        assert_eq!(allele_at(&read, 209), Some(b'T'), "block2 end -> seq[9]");
        assert_eq!(allele_at(&read, 210), None, "past read -> None");
    }

    /// "3M2D3M": a deletion. Positions inside the D return None; positions after map correctly
    /// (ref advanced over the deletion, read did NOT).
    #[test]
    fn bridge_allele_at_deletion() {
        // block1: ref 50..53 <-> seq[0..3] = "GGG"
        // del:    ref 53..55 (D2), consumes no read
        // block2: ref 55..58 <-> seq[3..6] = "TAC"
        let read = aligned(50, &[('M', 3), ('D', 2), ('M', 3)], b"GGGTAC");
        assert_eq!(allele_at(&read, 50), Some(b'G'));
        assert_eq!(allele_at(&read, 52), Some(b'G'), "last base before del -> seq[2]");
        assert_eq!(allele_at(&read, 53), None, "first deleted ref pos -> None");
        assert_eq!(allele_at(&read, 54), None, "second deleted ref pos -> None");
        assert_eq!(allele_at(&read, 55), Some(b'T'), "after del ref 55 -> seq[3]");
        assert_eq!(allele_at(&read, 57), Some(b'C'), "after del ref 57 -> seq[5]");
        assert_eq!(allele_at(&read, 58), None, "past read -> None");
    }

    /// "4S10M": softclip consumes read seq but not ref. The first matched ref position
    /// returns seq[4] (the 4 soft-clipped bases are skipped on the read axis).
    #[test]
    fn bridge_allele_at_softclip() {
        // softclip: seq[0..4] = "XXXX" (not aligned)
        // block:    ref 1000..1010 <-> seq[4..14] = "ACGTACGTAC"
        let read = aligned(1000, &[('S', 4), ('M', 10)], b"XXXXACGTACGTAC");
        assert_eq!(allele_at(&read, 1000), Some(b'A'), "first matched ref -> seq[4]");
        assert_eq!(allele_at(&read, 1001), Some(b'C'), "-> seq[5]");
        assert_eq!(allele_at(&read, 1009), Some(b'C'), "last matched -> seq[13]");
        assert_eq!(allele_at(&read, 999), None, "before ref_start -> None");
        assert_eq!(allele_at(&read, 1010), None, "past the end -> None");
    }

    /// "5M2I5M": an insertion. Ref positions are continuous across the insertion; the base
    /// after the insertion uses an advanced read coordinate (read += inserted bases).
    #[test]
    fn bridge_allele_at_insertion() {
        // block1: ref 10..15 <-> seq[0..5] = "AAAAA"
        // ins:    seq[5..7] = "II" (consumes read, no ref)
        // block2: ref 15..20 <-> seq[7..12] = "CCCCG"
        let read = aligned(10, &[('M', 5), ('I', 2), ('M', 5)], b"AAAAAIICCCCG");
        assert_eq!(allele_at(&read, 14), Some(b'A'), "last base of block1 -> seq[4]");
        // ref is continuous: 15 is the first base of block2, read coord jumped over the 2 ins.
        assert_eq!(allele_at(&read, 15), Some(b'C'), "after insertion ref 15 -> seq[7]");
        assert_eq!(allele_at(&read, 19), Some(b'G'), "block2 end ref 19 -> seq[11]");
        assert_eq!(allele_at(&read, 20), None, "past read -> None");
    }

    /// intron_chain_of("5M100N5M") == [(ref_start+5, ref_start+105)].
    #[test]
    fn bridge_intron_chain_single() {
        let read = aligned(100, &[('M', 5), ('N', 100), ('M', 5)], b"AAAAACCCCC");
        assert_eq!(intron_chain_of(&read), vec![(105, 205)]);
    }

    /// intron_chain over two introns, with a deletion that must NOT split the chain.
    #[test]
    fn bridge_intron_chain_multi_with_deletion() {
        // ref 100: 3M (100..103), 50N (103..153), 2M(153..155), 1D(155..156),
        //          2M(156..158), 20N(158..178), 4M(178..182)
        let read = aligned(
            100,
            &[
                ('M', 3),
                ('N', 50),
                ('M', 2),
                ('D', 1),
                ('M', 2),
                ('N', 20),
                ('M', 4),
            ],
            b"AAACCGGTTTT",
        );
        assert_eq!(
            intron_chain_of(&read),
            vec![(103, 153), (158, 178)],
            "two introns; the deletion does not break the chain"
        );
    }

    /// End-to-end build_read_obs: a spliced read + 3 PSV positions (exon1 / intron / exon2)
    /// -> ReadObs with the right intron_chain and psv_alleles = [Some, None, Some].
    #[test]
    fn bridge_build_read_obs_end_to_end() {
        // block1: ref 100..105 <-> seq[0..5] = "ACGTA"
        // intron: ref 105..205 (N100)
        // block2: ref 205..210 <-> seq[5..10] = "TTGCA"
        let read = aligned(100, &[('M', 5), ('N', 100), ('M', 5)], b"ACGTATTGCA");
        // PSV positions: one in exon1 (102 -> seq[2]='G'), one in the intron (150 -> None),
        // one in exon2 (207 -> seq[7]='G').
        let psv = [102u64, 150u64, 207u64];
        let obs = build_read_obs(&read, &psv);
        assert_eq!(obs.intron_chain, vec![(105, 205)], "intron chain from the N op");
        assert_eq!(
            obs.psv_alleles,
            vec![Some(b'G'), None, Some(b'G')],
            "[exon1 base, intron None, exon2 base]"
        );
    }
}
