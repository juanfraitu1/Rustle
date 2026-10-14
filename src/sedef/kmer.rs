//! Rolling canonical k-mer hashing over DNA sequences.
//!
//! Bases are 2-bit encoded (A=0, C=1, G=2, T=3). k-mers containing any
//! non-ACGT base are skipped. The canonical form of a k-mer is
//! min(fwd_hash, revcomp_hash), so a segment and its reverse complement
//! produce the same hashes -- inverted duplications are found for free.
//! Final hashes are avalanched with splitmix64 so sequence-similar k-mers
//! do not collide into similar hash values (needed for bottom-sketches).

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[inline]
fn encode_base(b: u8) -> Option<u64> {
    const LUT: [i8; 256] = {
        let mut t = [-1i8; 256];
        t[b'A' as usize] = 0;
        t[b'a' as usize] = 0;
        t[b'C' as usize] = 1;
        t[b'c' as usize] = 1;
        t[b'G' as usize] = 2;
        t[b'g' as usize] = 2;
        t[b'T' as usize] = 3;
        t[b't' as usize] = 3;
        t[b'U' as usize] = 3;
        t[b'u' as usize] = 3;
        t
    };
    let v = LUT[b as usize];
    if v < 0 {
        None
    } else {
        Some(v as u64)
    }
}

/// Calls `f(start, hash)` for every canonical k-mer in `seq`.
/// `start` is the 0-based start coordinate of the k-mer.
pub fn for_each_kmer(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64)) {
    for_each_kmer_raw(seq, k, |pos, fwd, rev| {
        f(pos, splitmix64(fwd.min(rev)));
    });
}

/// Like `for_each_kmer` but emits the RAW forward and revcomp hashes
/// (unmixed, not canonicalised). Stage 3 uses this for strand voting
/// (fwd-vs-rev hit counts) and anchor extraction.
pub fn for_each_kmer_raw(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64, u64)) {
    if k == 0 || k > 32 || seq.len() < k {
        return;
    }
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut run = 0usize; // consecutive valid bases
    for (i, &b) in seq.iter().enumerate() {
        match encode_base(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << (2 * k - 2));
                run += 1;
            }
            None => {
                run = 0;
                fwd = 0;
                rev = 0;
            }
        }
        if run >= k {
            f(i + 1 - k, fwd, rev);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_is_strand_symmetric() {
        let fwd = b"ACGTTGCA";
        let rev = b"TGCAACGT"; // reverse complement
        let h1: Vec<u64> = {
            let mut v = Vec::new();
            for_each_kmer(fwd, 4, |_, h| v.push(h));
            v
        };
        let h2: Vec<u64> = {
            let mut v = Vec::new();
            for_each_kmer(rev, 4, |_, h| v.push(h));
            v
        };
        let mut h1 = h1;
        let mut h2 = h2;
        h1.sort_unstable();
        h2.sort_unstable();
        assert_eq!(h1, h2);
    }

    #[test]
    fn skips_ambiguous_bases() {
        let mut n = 0;
        for_each_kmer(b"ACGTNNNNACGT", 4, |_, _| n += 1);
        assert_eq!(n, 2);
    }
}
