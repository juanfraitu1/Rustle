//! Augment-and-linearize certificate for reference-absent copies (v2.1). Pure statistics + a deterministic
//! dinucleotide-preserving decoy; the minimap2 re-alignment is injected as a closure so this module is testable
//! without any subprocess.

/// Deterministic LCG for reproducible shuffles (no external rng crate; Date/rand not needed).
struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self { Lcg(seed ^ 0x9E37_79B9_7F4A_7C15) }
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); self.0 }
    fn below(&mut self, n: usize) -> usize { if n <= 1 { 0 } else { (self.next() >> 33) as usize % n } }
}
fn fisher_yates<T>(v: &mut [T], rng: &mut Lcg) {
    for i in (1..v.len()).rev() { let j = rng.below(i + 1); v.swap(i, j); }
}

/// Altschul-Erikson (1985) dinucleotide-preserving shuffle via a random Eulerian path through the
/// nucleotide-transition multigraph. Same length, first/last base, and exact dinucleotide counts; scrambled
/// order. Deterministic in `seed`. Degenerate/short sequences (<3) or no valid Euler path -> a copy (a decoy
/// equal to the real candidate is conservative: it inflates the decoy fraction, never the real-minus-decoy gap).
pub fn dinucleotide_shuffle(seq: &[u8], seed: u64) -> Vec<u8> {
    use std::collections::BTreeMap;
    let n = seq.len();
    if n < 3 { return seq.to_vec(); }
    let last = seq[n - 1];
    let mut base_edges: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
    for w in seq.windows(2) { base_edges.entry(w[0]).or_default().push(w[1]); }
    let nodes: Vec<u8> = base_edges.keys().copied().collect();
    let mut rng = Lcg::new(seed);
    for _attempt in 0..1000 {
        // shuffle each node's outgoing edges; the LAST edge of each node (!= terminal) is its "tree" edge.
        let mut edges = base_edges.clone();
        for v in edges.values_mut() { fisher_yates(v, &mut rng); }
        let tree: BTreeMap<u8, u8> = edges.iter().filter(|(&b, _)| b != last)
            .map(|(&b, v)| (b, *v.last().unwrap())).collect();
        // valid iff following tree edges from every node reaches `last` without cycling (arborescence into last)
        let ok = nodes.iter().all(|&start| {
            if start == last { return true; }
            let (mut cur, mut steps) = (start, 0usize);
            loop {
                match tree.get(&cur) {
                    Some(&nx) => { cur = nx; if cur == last { break true; } steps += 1; if steps > nodes.len() { break false; } }
                    None => break false,
                }
            }
        });
        if !ok { continue; }
        // traverse: consume each node's edges in order (tree edge, at the end, is used last)
        let mut cursor: BTreeMap<u8, usize> = nodes.iter().map(|&b| (b, 0)).collect();
        let mut out = Vec::with_capacity(n);
        out.push(seq[0]);
        let mut cur = seq[0];
        while out.len() < n {
            let v = match edges.get(&cur) { Some(v) => v, None => break };
            let i = cursor[&cur];
            if i >= v.len() { break; }
            let nx = v[i];
            *cursor.get_mut(&cur).unwrap() += 1;
            out.push(nx);
            cur = nx;
        }
        if out.len() == n { return out; }
    }
    seq.to_vec()
}

/// Verdict from linearize_certificate: the candidate linearizes relative to null/decoy shuffles.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Verdict {
    Linearizes,
    Not,
    Undetermined,
}

/// Result of the linearize_certificate test: candidate fraction, decoy statistics, and verdict.
#[derive(Clone, Debug)]
pub struct LinearizeCertificate {
    pub n_pool: usize,
    pub linearized_frac_real: f64,
    pub mean_frac_decoy: f64,
    pub delta: f64,
    pub perm_p: f64,
    pub verdict: Verdict,
}

/// Fraction of reads whose primary hit is the candidate contig with MAPQ > 0.
fn frac_on_candidate(hits: &[Option<(usize, u32)>], cand_idx: usize) -> f64 {
    if hits.is_empty() {
        return 0.0;
    }
    let k = hits
        .iter()
        .filter(|h| matches!(h, Some((i, q)) if *i == cand_idx && *q > 0))
        .count();
    k as f64 / hits.len() as f64
}

/// Test whether a sequence (candidate) linearizes: its primary-with-MAPQ>0 fraction significantly
/// exceeds the mean of N dinucleotide-shuffled decoys + a reverse-complement decoy. Pure function;
/// the realign closure is injected for testability.
///
/// # Arguments
/// - `candidate_seq`: the candidate contig (appended as the last element to family_copy_seqs).
/// - `family_copy_seqs`: the background copies (NOT including the candidate).
/// - `pool_reads`: the read pool for evaluation.
/// - `n_decoys`: count of dinucleotide-shuffled decoys.
/// - `seed`: RNG seed for deterministic shuffles.
/// - `min_pool`: minimum pool size; if `n_pool < min_pool`, returns Undetermined with NaN fields.
/// - `alpha`: permutation p-value threshold (e.g., 0.05) for the Linearizes verdict.
/// - `realign`: injected closure `Fn(refs, reads) -> Vec<Option<(contig_idx, mapq)>>` per read.
pub fn linearize_certificate(
    candidate_seq: &[u8],
    family_copy_seqs: &[Vec<u8>],
    pool_reads: &[Vec<u8>],
    n_decoys: usize,
    seed: u64,
    min_pool: usize,
    alpha: f64,
    realign: impl Fn(&[Vec<u8>], &[Vec<u8>]) -> Vec<Option<(usize, u32)>>,
) -> LinearizeCertificate {
    let n_pool = pool_reads.len();
    if n_pool < min_pool {
        return LinearizeCertificate {
            n_pool,
            linearized_frac_real: f64::NAN,
            mean_frac_decoy: f64::NAN,
            delta: f64::NAN,
            perm_p: f64::NAN,
            verdict: Verdict::Undetermined,
        };
    }

    let cand_idx = family_copy_seqs.len(); // candidate is the LAST contig
    let build = |extra: &[u8]| -> Vec<Vec<u8>> {
        let mut v = family_copy_seqs.to_vec();
        v.push(extra.to_vec());
        v
    };

    // Compute the real candidate's linearized fraction.
    let real = frac_on_candidate(&realign(&build(candidate_seq), pool_reads), cand_idx);

    // Generate decoys: N dinucleotide shuffles (distinct seeds) + reverse-complement.
    let mut decoy_fracs: Vec<f64> = Vec::with_capacity(n_decoys + 1);
    for d in 0..n_decoys {
        let decoy = dinucleotide_shuffle(candidate_seq, seed.wrapping_add(d as u64 + 1));
        decoy_fracs.push(frac_on_candidate(&realign(&build(&decoy), pool_reads), cand_idx));
    }
    let rc = crate::vg_family::seq_utils::reverse_complement(candidate_seq);
    decoy_fracs.push(frac_on_candidate(
        &realign(&build(&rc), pool_reads),
        cand_idx,
    ));

    // Compute statistics.
    let nd = decoy_fracs.len();
    let mean_decoy = decoy_fracs.iter().sum::<f64>() / nd as f64;
    let n_ge = decoy_fracs.iter().filter(|&&d| d >= real).count();
    let perm_p = (n_ge as f64 + 1.0) / (nd as f64 + 1.0);
    let verdict = if perm_p < alpha {
        Verdict::Linearizes
    } else {
        Verdict::Not
    };

    LinearizeCertificate {
        n_pool,
        linearized_frac_real: real,
        mean_frac_decoy: mean_decoy,
        delta: real - mean_decoy,
        perm_p,
        verdict,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn dinuc_counts(s: &[u8]) -> std::collections::BTreeMap<(u8,u8), usize> {
        let mut m = std::collections::BTreeMap::new();
        for w in s.windows(2) { *m.entry((w[0], w[1])).or_insert(0) += 1; }
        m
    }
    #[test]
    fn dinucleotide_shuffle_preserves_composition_and_is_deterministic() {
        let seq = b"ACGTACGTTTGGCCAAACGTACGTGGGCCCAAATTT";
        let a = dinucleotide_shuffle(seq, 42);
        let b = dinucleotide_shuffle(seq, 42);
        assert_eq!(a, b, "deterministic for a given seed");
        assert_eq!(a.len(), seq.len(), "length preserved");
        assert_eq!(a[0], seq[0], "first base preserved");
        assert_eq!(*a.last().unwrap(), *seq.last().unwrap(), "last base preserved");
        assert_eq!(dinuc_counts(&a), dinuc_counts(seq), "exact dinucleotide counts preserved");
        let c = dinucleotide_shuffle(seq, 43);
        assert_ne!(a, c, "different seed -> different shuffle (for a shufflable seq)");
    }

    // A fake realign: a read "belongs" to the candidate iff its bytes equal the candidate contig's bytes;
    // then it maps uniquely (mapq 60) to the candidate index; otherwise it maps to copy 0 with mapq 0 (tied).
    fn fake_realign(refs: &[Vec<u8>], reads: &[Vec<u8>]) -> Vec<Option<(usize, u32)>> {
        let cand_idx = refs.len() - 1;
        reads
            .iter()
            .map(|r| {
                if r == &refs[cand_idx] {
                    Some((cand_idx, 60))
                } else {
                    Some((0, 0))
                }
            })
            .collect()
    }

    #[test]
    fn real_copy_linearizes_decoy_does_not() {
        // Use a longer sequence to minimize chance of shuffle returning the original
        let cand = b"ACGTACGTTTGGCCAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let copies = vec![b"TTTTGGGGCCCCAAAAAAAATTTTGGGG".to_vec()];
        // 8 pool reads that ARE the candidate (they linearize on it), 2 that are noise
        let mut pool: Vec<Vec<u8>> = (0..8).map(|_| cand.clone()).collect();
        pool.push(b"NNNNNN".to_vec());
        pool.push(b"NNNNNN".to_vec());
        let cert = linearize_certificate(&cand, &copies, &pool, 19, 7, 5, 0.05, fake_realign);
        assert!(
            (cert.linearized_frac_real - 0.8).abs() < 1e-9,
            "8/10 land on candidate"
        );
        assert!(cert.mean_frac_decoy == 0.0, "decoys != candidate bytes -> no read lands on them");
        assert!(cert.delta > 0.5);
        assert!(
            cert.perm_p <= 1.0 / 20.0 + 1e-9,
            "no decoy beats real -> perm_p = 1/(N+1)"
        );
        assert!(matches!(cert.verdict, Verdict::Linearizes));
    }

    #[test]
    fn null_candidate_is_not_linearized() {
        let cand = b"ACGTACGTTTGGCCAAACGTACGT".to_vec();
        let copies = vec![b"TTTTGGGGCCCCAAAA".to_vec()];
        let pool: Vec<Vec<u8>> = (0..10).map(|_| b"NNNNNN".to_vec()).collect(); // nothing matches candidate
        let cert = linearize_certificate(&cand, &copies, &pool, 19, 7, 5, 0.05, fake_realign);
        assert_eq!(cert.linearized_frac_real, 0.0);
        assert!(cert.perm_p > 0.05, "real == decoys (both 0) -> perm_p large");
        assert!(matches!(cert.verdict, Verdict::Not));
    }

    #[test]
    fn small_pool_is_undetermined() {
        let cand = b"ACGT".to_vec();
        let pool = vec![cand.clone(), cand.clone()];
        let cert = linearize_certificate(&cand, &[], &pool, 19, 7, 5, 0.05, fake_realign);
        assert!(matches!(cert.verdict, Verdict::Undetermined));
    }
}
