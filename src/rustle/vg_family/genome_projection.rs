//! Liftoff-style genome-projection copy enumeration (spec §7): project a de-novo family consensus onto the
//! genome to enumerate near-identical genomic copies (famCN), recovering K=0 collapses RNA merges. In-engine
//! minimap2 (no Liftoff dependency); seeded by our own consensus, so no reference-annotation circularity.
use anyhow::Result;
use std::io::Write;

#[derive(Clone, Debug)]
pub struct CopyLocus { pub chrom: String, pub start: u64, pub end: u64, pub identity: f64 }

/// minimap2 the consensus against the genome; keep hits with identity ≥ `min_identity`, aligned-fraction
/// of the consensus ≥ `min_cov` (structure-preserving), disjoint from each other and from `known` loci.
pub fn project_family_copies(
    consensus: &[u8],
    genome_fasta: &str,
    known: &[(String, u64, u64)],
    min_identity: f64,
    min_cov: f64,
    minimap2: &str,
    threads: usize,
) -> Result<Vec<CopyLocus>> {
    let dir = std::env::temp_dir();
    let q = dir.join(format!("rustle_proj_q_{}_{}.fa", std::process::id(), consensus.len()));
    struct Cl(std::path::PathBuf); impl Drop for Cl { fn drop(&mut self) { let _ = std::fs::remove_file(&self.0); } }
    let _c = Cl(q.clone());
    { let mut f = std::fs::File::create(&q)?; writeln!(f, ">cons")?; f.write_all(consensus)?; writeln!(f)?; }
    let out = std::process::Command::new(minimap2)
        // -p 0.01: report divergent secondaries too (default -x splice -p suppresses them, hiding all but
        // near-identical copies); the id/cov filter below decides which to keep.
        .args(["-c", "-x", "splice", "-N", "50", "-p", "0.01", "-t"]).arg(threads.to_string())
        .arg(genome_fasta).arg(&q).output()
        .map_err(|e| anyhow::anyhow!("minimap2 ('{minimap2}') projection failed: {e}"))?;
    if !out.status.success() { return Ok(Vec::new()); }
    let clen = consensus.len() as f64;
    let mut hits: Vec<CopyLocus> = Vec::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let tname = f[5].to_string();
        let ts = f[7].parse::<u64>().unwrap_or(0);
        let te = f[8].parse::<u64>().unwrap_or(0);
        let qs = f[2].parse::<f64>().unwrap_or(0.0);
        let qe = f[3].parse::<f64>().unwrap_or(0.0);
        let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = de.map(|d| 1.0 - d).unwrap_or_else(|| {
            f[9].parse::<f64>().unwrap_or(0.0) / f[10].parse::<f64>().unwrap_or(1.0).max(1.0)
        });
        let cov = (qe - qs) / clen.max(1.0);
        if ident >= min_identity && cov >= min_cov {
            hits.push(CopyLocus { chrom: tname, start: ts, end: te, identity: ident });
        }
    }
    // disjoint filter: sort, drop hits overlapping an already-kept hit or a `known` locus.
    hits.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
    let mut kept: Vec<CopyLocus> = Vec::new();
    let overlaps = |c: &CopyLocus, k: &(String, u64, u64)| c.chrom == k.0 && c.start < k.2 && k.1 < c.end;
    for h in hits {
        if kept.iter().any(|k| k.chrom == h.chrom && k.start < h.end && h.start < k.end) { continue; }
        if known.iter().any(|k| overlaps(&h, k)) { continue; }
        kept.push(h);
    }
    Ok(kept)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projection_enumerates_disjoint_copies() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        // Build a tiny 2-contig genome file with the SAME 1kb sequence at two loci = 2 genomic copies.
        let dir = std::env::temp_dir().join(format!("rustle_proj_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // NOTE: a naive affine index formula like `(i*7+3)%4` is periodic with period 4 (any `i*k+c mod 4`
        // with k odd cycles every 4 bases) -- i.e. literally "TGCA" repeated. minimap2 then finds a match at
        // every 4bp-phase-aligned offset across the whole c1 contig, bridging the two real tandem copies into
        // one contiguous disjoint-kept interval (verified: 51 raw PAF hits -> 2 kept, not 3). Use a splitmix64
        // hash per index instead so the sequence is non-degenerate, deterministic, and dependency-free.
        let splitmix = |i: u64| -> u64 {
            let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        };
        let seq: String = (0..1000u64).map(|i| "ACGT".as_bytes()[(splitmix(i) % 4) as usize] as char).collect();
        let fa = dir.join("g.fa");
        std::fs::write(&fa, format!(">c1\n{seq}{seq}\n>c2\n{seq}\n")).unwrap(); // c1 has 2 tandem copies, c2 has 1
        let hits = project_family_copies(seq.as_bytes(), fa.to_str().unwrap(), &[], 0.98, 0.90, "minimap2", 2).unwrap();
        assert!(hits.len() >= 3, "3 genomic copies (2 on c1, 1 on c2) expected, got {}", hits.len());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Regression for the `-p` (secondary-to-primary score ratio) bug: `minimap2 -x splice`'s DEFAULT `-p`
    /// suppresses divergent secondary alignments whose chain overlaps the (higher-scoring) near-identical
    /// hit in query space -- so a family with one 0%-divergent copy and two divergent copies (~8%, ~15%)
    /// projects to only the identical copy, silently dropping the divergent ones. `-p 0.01` makes minimap2
    /// report them; the existing identity/coverage filter then decides which to keep.
    #[test]
    fn projection_finds_divergent_copies_not_just_identical() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }

        // Deterministic splitmix64-based non-degenerate sequence generator (see note on the test above:
        // avoid periodic index formulas that alias into a repeated motif and confuse minimap2).
        let splitmix = |i: u64| -> u64 {
            let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        };
        let bases = [b'A', b'C', b'G', b'T'];
        let gen_seq = |seed: u64, len: u64| -> Vec<u8> {
            (0..len)
                .map(|i| bases[(splitmix(seed.wrapping_mul(0x2545_F491_4F6C_DD1D).wrapping_add(i)) % 4) as usize])
                .collect::<Vec<u8>>()
        };
        // Mutate an EXACT `frac` of positions to a different base (rank-selected by hash, not a per-base
        // coin flip) so the realized divergence matches `frac` precisely regardless of sequence length.
        let mutate = |seq: &[u8], frac: f64, seed: u64| -> Vec<u8> {
            let n = seq.len();
            let mut idx: Vec<usize> = (0..n).collect();
            idx.sort_by_key(|&i| splitmix(seed.wrapping_add(i as u64)));
            let n_mut = (frac * n as f64).round() as usize;
            let mut out = seq.to_vec();
            for &i in idx.iter().take(n_mut) {
                let orig = out[i];
                let h = splitmix(seed ^ 0xABCDEF ^ i as u64);
                let mut alt = bases[(h % 4) as usize];
                if alt == orig { alt = bases[((h % 4) as usize + 1) % 4]; }
                out[i] = alt;
            }
            out
        };

        let copy_len = 1500u64;
        let copy0 = gen_seq(1, copy_len); // 0% divergent: the consensus itself
        let copy8 = mutate(&copy0, 0.08, 100); // ~8% divergent
        let copy15 = mutate(&copy0, 0.15, 200); // ~15% divergent

        // Assemble one contig: bg + copy0 + ~20kb bg + copy8 + ~20kb bg + copy15 + bg, so the three copies
        // sit at distinct loci ~20kb apart, not adjacent/tandem.
        let dir = std::env::temp_dir().join(format!("rustle_proj_div_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut genome: Vec<u8> = Vec::new();
        genome.extend_from_slice(&gen_seq(9001, 2_000));
        genome.extend_from_slice(&copy0);
        genome.extend_from_slice(&gen_seq(9002, 20_000));
        genome.extend_from_slice(&copy8);
        genome.extend_from_slice(&gen_seq(9003, 20_000));
        genome.extend_from_slice(&copy15);
        genome.extend_from_slice(&gen_seq(9004, 2_000));

        let fa = dir.join("g.fa");
        std::fs::write(&fa, format!(">chr1\n{}\n", String::from_utf8(genome).unwrap())).unwrap();

        let hits = project_family_copies(&copy0, fa.to_str().unwrap(), &[], 0.80, 0.50, "minimap2", 2).unwrap();
        assert!(
            hits.len() >= 3,
            "expected >=3 disjoint copies (identical + ~8% + ~15% divergent), got {}: {:?}",
            hits.len(),
            hits
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
