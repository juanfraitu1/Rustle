//! Small sequence utilities the family-analysis modules depend on.
//!
//! Two reverse complements with DIFFERENT semantics live here, named by what they do:
//! - `reverse_complement` — uppercase ACGT only; every other byte (lowercase included) maps to `N`.
//!   Relocated verbatim from the retired assembler `vg.rs` (`docs/RETIREMENT_AND_MIGRATION.md`).
//! - `revcomp_keep_case` — complements both cases (`N`/`n` kept), passes every other byte through
//!   unchanged (Python `str.translate` semantics).
//!
//! plus `hw_distance` / `aln_id`, the HW (infix) edit distance and identity that equal edlib's
//! `mode="HW"` bit for bit. `revcomp_keep_case`, `hw_distance` and `aln_id` were moved here from
//! `bridge_detector.rs` (a port of `bench/recombination_bridge_detector.py`) when the rest of that module
//! was removed as dead code (2026-09-24, tag `notebook-2026-09-24`).
//!
//! **STATUS:** INFRASTRUCTURE  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

/// Reverse-complement a nucleotide sequence. Non-ACGT bytes map to `b'N'`.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// One-base complement, exactly the Python `_COMP = str.maketrans("ACGTacgtNn",
/// "TGCAtgcaNn")`: A<->T, C<->G, G<->C, T<->A (both cases), N->N, n->n; EVERY OTHER
/// byte is left unchanged (Python `str.translate` passes unmapped chars through).
#[inline]
fn comp_base_keep_case(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

/// Case-preserving reverse complement: complement every base via `comp_base_keep_case`, then reverse
/// (`s.translate(_COMP)[::-1]`).
pub fn revcomp_keep_case(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp_base_keep_case(b)).collect()
}

/// HW (infix) edit distance of `q` against `t`: the minimum edit distance of `q` to ANY
/// substring of `t` (free gaps at both ends of `t`). Equals `edlib.align(q, t,
/// mode="HW", task="distance")["editDistance"]`. Two rolling DP rows over `t` columns
/// (`dp[0][*]=0`, `dp[i][0]=i`, answer `= min_j dp[|q|][j]`).
pub(crate) fn hw_distance(q: &[u8], t: &[u8]) -> usize {
    let lq = q.len();
    let lt = t.len();
    // dp row 0 = 0 across all t columns (a free start position in t).
    let mut prev: Vec<usize> = vec![0; lt + 1];
    let mut cur: Vec<usize> = vec![0; lt + 1];
    for i in 1..=lq {
        cur[0] = i; // dp[i][0] = i
        let qi = q[i - 1];
        for j in 1..=lt {
            let sub = prev[j - 1] + if qi == t[j - 1] { 0 } else { 1 };
            let del = prev[j] + 1;
            let ins = cur[j - 1] + 1;
            cur[j] = sub.min(del).min(ins);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    // answer = min over j of dp[|q|][j]  (free end position in t). After the final swap
    // the last computed row is in `prev`. For lq == 0, prev is the all-zero row 0.
    *prev.iter().min().expect("row has >= 1 column")
}

/// Best HW (infix) identity of `q` inside `t`, trying `q` AND `revcomp_keep_case(q)`.
/// `id = 1 - min_dist / len(q)`; `0.0` if either string is empty
/// (`recombination_bridge_detector.py:70`).
pub fn aln_id(q: &[u8], t: &[u8]) -> f64 {
    if q.is_empty() || t.is_empty() {
        return 0.0;
    }
    let d_fwd = hw_distance(q, t);
    let d_rev = hw_distance(&revcomp_keep_case(q), t);
    let best = d_fwd.min(d_rev);
    // 1.0 - best / len(q)  (single f64 division, matches Python exactly).
    1.0 - best as f64 / q.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn reverse_complement_matches_legacy_semantics() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
        // non-ACGT (incl. lowercase) -> N; reverse THEN complement: [A,N,a] -> [a,N,A] -> [N,N,T]
        assert_eq!(reverse_complement(b"ANa"), b"NNT");
    }

    #[test]
    fn revcomp_table_semantics() {
        // A<->T, C<->G, N->N (both cases), everything else passes through, then reverse.
        assert_eq!(revcomp_keep_case(b"ACGT"), b"ACGT".to_vec()); // palindrome
        assert_eq!(revcomp_keep_case(b"AAAA"), b"TTTT".to_vec());
        assert_eq!(revcomp_keep_case(b"ACGTN"), b"NACGT".to_vec()); // N->N, reversed
        assert_eq!(revcomp_keep_case(b"acgt"), b"acgt".to_vec()); // lowercase palindrome
        assert_eq!(revcomp_keep_case(b"acg"), b"cgt".to_vec()); // a->t,c->g,g->c reversed
        assert_eq!(revcomp_keep_case(b""), Vec::<u8>::new());
        // unmapped byte passes through unchanged (Python str.translate leaves it)
        assert_eq!(revcomp_keep_case(b"AXT"), b"AXT".to_vec()); // T->A, X->X, A->T ; reversed = A X T
    }

    #[test]
    fn hw_distance_corners() {
        assert_eq!(hw_distance(b"ACGT", b"TTACGTGG"), 0); // exact infix
        assert_eq!(hw_distance(b"ACGT", b"ACGT"), 0);
        assert_eq!(hw_distance(b"ACGTACGTACGT", b"ACGT"), 8); // q longer -> 8 deletions
        assert_eq!(hw_distance(b"ACGT", b"ACAT"), 1); // one substitution
        assert_eq!(hw_distance(b"A", b"TTTT"), 1); // best infix "" or one mismatch
    }

    /// PROVES the pure-Rust HW-DP `aln_id` == edlib on >= 80 adversarial pairs (exact
    /// bit-for-bit float equality). Reports the first mismatch. The fixture is the `aln_id` section of
    /// the retired `bridge_detector_fixture.json`.
    #[test]
    fn aln_id_parity_vs_edlib() {
        let fx: serde_json::Value = serde_json::from_str(include_str!("testdata/aln_id_fixture.json"))
            .expect("parse aln_id fixture json");
        let cases = fx["aln_id"].as_array().expect("aln_id array");
        assert!(cases.len() >= 80, "need >= 80 aln_id cases, got {}", cases.len());
        let mut classes: BTreeSet<String> = BTreeSet::new();
        let mut n_empty = 0;
        let mut n_n = 0;
        let mut n_revcomp = 0;
        for c in cases {
            let q = c["q"].as_str().unwrap();
            let t = c["t"].as_str().unwrap();
            // EXACT bits (serde_json's decimal float parser is 1-ULP imprecise).
            let want = f64::from_bits(c["bits"].as_u64().unwrap());
            let cls = c["cls"].as_str().unwrap_or("");
            classes.insert(cls.to_string());
            if q.is_empty() || t.is_empty() {
                n_empty += 1;
            }
            if q.contains('N') || t.contains('N') {
                n_n += 1;
            }
            if cls.contains("revcomp") {
                n_revcomp += 1;
            }
            let got = aln_id(q.as_bytes(), t.as_bytes());
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "aln_id MISMATCH cls='{cls}' q='{q}' t='{t}': rust {got} (bits {:#x}) != python {want} (bits {:#x})",
                got.to_bits(),
                want.to_bits()
            );
        }
        // edge-class coverage guard
        for required in ["both_empty", "q_empty", "t_empty", "q_longer", "exact_infix", "tandem"] {
            assert!(classes.contains(required), "aln_id fixture missing class '{required}'");
        }
        assert!(n_empty >= 3 && n_n >= 3 && n_revcomp >= 3, "aln_id fixture lacks edge coverage");
    }
}
