//! Small sequence utilities the family-analysis modules depend on.
//!
//! `reverse_complement` was relocated here verbatim from the retired assembler `vg.rs` during the
//! StringTie-assembler carve (`docs/RETIREMENT_AND_MIGRATION.md`): it is the one general symbol the
//! thesis layer imported from the assembler at runtime, so keeping it in a thesis-owned module lets the
//! whole assembler island be deleted. Behaviour is unchanged (non-ACGT bytes map to `N`, uppercase only).
//!
//! **STATUS:** INFRASTRUCTURE  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

/// Reverse-complement a nucleotide sequence. Non-ACGT bytes map to `b'N'`.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_matches_legacy_semantics() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
        // non-ACGT (incl. lowercase) -> N; reverse THEN complement: [A,N,a] -> [a,N,A] -> [N,N,T]
        assert_eq!(reverse_complement(b"ANa"), b"NNT");
    }
}
