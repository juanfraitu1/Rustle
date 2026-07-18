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
}
