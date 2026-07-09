//! Liftoff-style genome-projection copy enumeration (spec §7): project a de-novo family consensus onto the
//! genome to enumerate near-identical genomic copies (famCN), recovering K=0 collapses RNA merges. In-engine
//! minimap2 (no Liftoff dependency); seeded by our own consensus, so no reference-annotation circularity.
use anyhow::Result;
use std::collections::HashMap;
use std::io::Write;

#[derive(Clone, Debug)]
pub struct CopyLocus { pub chrom: String, pub start: u64, pub end: u64, pub identity: f64, pub cov: f64 }

/// RAII temp-file cleanup helper shared by the single and batch projection entry points.
struct TempFile(std::path::PathBuf);
impl Drop for TempFile { fn drop(&mut self) { let _ = std::fs::remove_file(&self.0); } }

/// Run `minimap2 -c -x splice -N 50 -p 0.01` (query FASTA at `query_path`) against `genome_fasta` and
/// return the raw PAF stdout. `-p 0.01`: report divergent secondaries too (default -x splice -p suppresses
/// them, hiding all but near-identical copies); the id/cov filter downstream decides which to keep. Returns
/// `Ok(None)` (not an error) if minimap2 exits non-zero, matching the existing graceful-degradation contract.
fn run_minimap2_paf(
    query_path: &std::path::Path,
    genome_fasta: &str,
    minimap2: &str,
    threads: usize,
) -> Result<Option<String>> {
    let out = std::process::Command::new(minimap2)
        .args(["-c", "-x", "splice", "-N", "50", "-p", "0.01", "-t"]).arg(threads.to_string())
        .arg(genome_fasta).arg(query_path).output()
        .map_err(|e| anyhow::anyhow!("minimap2 ('{minimap2}') projection failed: {e}"))?;
    if !out.status.success() { return Ok(None); }
    Ok(Some(String::from_utf8_lossy(&out.stdout).into_owned()))
}

/// Parse PAF text into per-query-name hit lists, filtered by identity/coverage. `qlen` (PAF field 1) is
/// read directly per-record, so this works uniformly whether the PAF came from a single-query or a
/// multi-query (batch) minimap2 run — no external query-length bookkeeping needed.
fn parse_paf_hits(paf: &str, min_identity: f64, min_cov: f64) -> HashMap<String, Vec<CopyLocus>> {
    let mut by_query: HashMap<String, Vec<CopyLocus>> = HashMap::new();
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let qname = f[0].to_string();
        let qlen = f[1].parse::<f64>().unwrap_or(1.0).max(1.0);
        let tname = f[5].to_string();
        let ts = f[7].parse::<u64>().unwrap_or(0);
        let te = f[8].parse::<u64>().unwrap_or(0);
        let qs = f[2].parse::<f64>().unwrap_or(0.0);
        let qe = f[3].parse::<f64>().unwrap_or(0.0);
        let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = de.map(|d| 1.0 - d).unwrap_or_else(|| {
            f[9].parse::<f64>().unwrap_or(0.0) / f[10].parse::<f64>().unwrap_or(1.0).max(1.0)
        });
        let cov = (qe - qs) / qlen;
        if ident >= min_identity && cov >= min_cov {
            by_query.entry(qname).or_default().push(CopyLocus { chrom: tname, start: ts, end: te, identity: ident, cov });
        }
    }
    by_query
}

/// Disjoint filter: sort by (chrom, start), drop hits overlapping an already-kept hit or a `known` locus.
fn disjoint_filter(mut hits: Vec<CopyLocus>, known: &[(String, u64, u64)]) -> Vec<CopyLocus> {
    hits.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
    let mut kept: Vec<CopyLocus> = Vec::new();
    let overlaps = |c: &CopyLocus, k: &(String, u64, u64)| c.chrom == k.0 && c.start < k.2 && k.1 < c.end;
    for h in hits {
        if kept.iter().any(|k| k.chrom == h.chrom && k.start < h.end && h.start < k.end) { continue; }
        if known.iter().any(|k| overlaps(&h, k)) { continue; }
        kept.push(h);
    }
    kept
}

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
    let _c = TempFile(q.clone());
    { let mut f = std::fs::File::create(&q)?; writeln!(f, ">cons")?; f.write_all(consensus)?; writeln!(f)?; }
    let paf = match run_minimap2_paf(&q, genome_fasta, minimap2, threads)? {
        Some(p) => p,
        None => return Ok(Vec::new()),
    };
    let mut by_query = parse_paf_hits(&paf, min_identity, min_cov);
    let hits = by_query.remove("cons").unwrap_or_default();
    Ok(disjoint_filter(hits, known))
}

/// Batch variant: ONE minimap2 invocation (one genome index load) for MANY families' consensuses, instead
/// of one invocation per family. Writes all consensuses to a single multi-record query FASTA (header =
/// `family_id`), runs minimap2 once, groups PAF hits by query name, then applies the same per-family
/// identity/coverage filter and disjoint+known-exclusion filter as `project_family_copies`.
pub fn project_families_batch(
    consensuses: &[(String, Vec<u8>)],
    genome_fasta: &str,
    known: &HashMap<String, Vec<(String, u64, u64)>>,
    min_identity: f64,
    min_cov: f64,
    minimap2: &str,
    threads: usize,
) -> Result<HashMap<String, Vec<CopyLocus>>> {
    let dir = std::env::temp_dir();
    let q = dir.join(format!("rustle_proj_batch_{}_{}.fa", std::process::id(), consensuses.len()));
    let _c = TempFile(q.clone());
    {
        let mut f = std::fs::File::create(&q)?;
        for (fam_id, seq) in consensuses {
            if seq.is_empty() { continue; }
            writeln!(f, ">{fam_id}")?;
            f.write_all(seq)?;
            writeln!(f)?;
        }
    }
    let empty: HashMap<String, Vec<CopyLocus>> = consensuses.iter().map(|(id, _)| (id.clone(), Vec::new())).collect();
    let paf = match run_minimap2_paf(&q, genome_fasta, minimap2, threads)? {
        Some(p) => p,
        None => return Ok(empty),
    };
    let by_query = parse_paf_hits(&paf, min_identity, min_cov);
    let no_known: Vec<(String, u64, u64)> = Vec::new();
    let mut result: HashMap<String, Vec<CopyLocus>> = HashMap::new();
    for (fam_id, _) in consensuses {
        let hits = by_query.get(fam_id).cloned().unwrap_or_default();
        let k = known.get(fam_id).map(|v| v.as_slice()).unwrap_or(no_known.as_slice());
        result.insert(fam_id.clone(), disjoint_filter(hits, k));
    }
    Ok(result)
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

    /// TDD for `project_families_batch`: ONE minimap2 pass must correctly enumerate + bucket-by-identity
    /// loci for MULTIPLE families sharing one query FASTA, without cross-family leakage.
    #[test]
    fn project_families_batch_one_pass_buckets_by_identity() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }

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
        // Family F1: sequence "A" at 0% / ~8% / ~15% divergence.
        let a0 = gen_seq(1, copy_len);
        let a8 = mutate(&a0, 0.08, 100);
        let a15 = mutate(&a0, 0.15, 200);
        // A SHORT ~50%-length near-identical fragment of A (first half, ~1% mutated): pins the fragment
        // concern -- it clears the totalCN floor (cov>=0.50) but must be EXCLUDED from famCN (cov<0.90).
        let a_frag = mutate(&a0[..(copy_len as usize / 2)], 0.01, 400);
        // Family F2: a DIFFERENT sequence "B" (distinct seed) at 0% / ~2% divergence.
        let b0 = gen_seq(500_001, copy_len);
        let b2 = mutate(&b0, 0.02, 300);

        let dir = std::env::temp_dir().join(format!("rustle_proj_batch_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut genome: Vec<u8> = Vec::new();
        genome.extend_from_slice(&gen_seq(9001, 2_000));
        genome.extend_from_slice(&a0);
        genome.extend_from_slice(&gen_seq(9002, 20_000));
        genome.extend_from_slice(&a8);
        genome.extend_from_slice(&gen_seq(9003, 20_000));
        genome.extend_from_slice(&a15);
        genome.extend_from_slice(&gen_seq(9007, 20_000));
        genome.extend_from_slice(&a_frag);
        genome.extend_from_slice(&gen_seq(9004, 20_000));
        genome.extend_from_slice(&b0);
        genome.extend_from_slice(&gen_seq(9005, 20_000));
        genome.extend_from_slice(&b2);
        genome.extend_from_slice(&gen_seq(9006, 2_000));

        let fa = dir.join("g.fa");
        std::fs::write(&fa, format!(">chr1\n{}\n", String::from_utf8(genome).unwrap())).unwrap();

        let consensuses: Vec<(String, Vec<u8>)> = vec![("F1".to_string(), a0.clone()), ("F2".to_string(), b0.clone())];
        let known: std::collections::HashMap<String, Vec<(String, u64, u64)>> = std::collections::HashMap::new();
        let result = project_families_batch(&consensuses, fa.to_str().unwrap(), &known, 0.80, 0.50, "minimap2", 2).unwrap();

        let f1 = result.get("F1").expect("F1 must be present in batch result");
        let f2 = result.get("F2").expect("F2 must be present in batch result");

        assert!(f1.len() >= 3, "F1 expected >=3 loci (0%/~8%/~15% divergent), got {}: {:?}", f1.len(), f1);
        assert!(f2.len() >= 2, "F2 expected >=2 loci (0%/~2% divergent), got {}: {:?}", f2.len(), f2);

        // Identity spans: F1 should carry a near-1.0, a ~0.92, and a ~0.85 hit.
        let f1_has = |lo: f64, hi: f64| f1.iter().any(|c| c.identity >= lo && c.identity <= hi);
        assert!(f1_has(0.97, 1.0), "F1 missing ~1.0-identity hit: {:?}", f1);
        assert!(f1_has(0.88, 0.96), "F1 missing ~0.92-identity hit: {:?}", f1);
        assert!(f1_has(0.80, 0.88), "F1 missing ~0.85-identity hit: {:?}", f1);

        // DISJOINT identity bands for F2, so the two asserts pin two distinct hits (not one): the
        // identical copy sits in [0.99,1.0], the ~2%-divergent copy in [0.96,0.99).
        let f2_has = |lo: f64, hi: f64| f2.iter().any(|c| c.identity >= lo && c.identity < hi);
        assert!(f2_has(0.99, 1.0001), "F2 missing ~1.0-identity hit: {:?}", f2);
        assert!(f2_has(0.96, 0.99), "F2 missing ~0.98-identity hit: {:?}", f2);

        // CopyLocus.cov must be populated: F1's full-length copies align (nearly) the whole query.
        let f1_full: Vec<&CopyLocus> = f1.iter().filter(|c| c.cov >= 0.90).collect();
        assert!(f1_full.len() >= 3,
            "F1 expected >=3 FULL-LENGTH copies (cov>=0.90), got {}: {:?}", f1_full.len(), f1);
        assert!(f1.iter().any(|c| c.identity >= 0.97 && c.cov >= 0.95),
            "F1's identical copy should have cov close to 1.0: {:?}", f1);

        // The planted ~50% fragment (cov~0.5, id~0.99) must appear (clears totalCN floor cov>=0.50) but
        // must be EXCLUDED from the famCN bucket (cov<0.90), even though its identity is >=0.98.
        let f1_famcn = f1.iter().filter(|c| c.identity >= 0.98 && c.cov >= 0.90).count();
        let f1_frag_hits = f1.iter().filter(|c| c.cov >= 0.40 && c.cov < 0.75).count();
        assert!(f1_frag_hits >= 1, "F1 half-length fragment (cov~0.5) should be present: {:?}", f1);
        // famCN counts only the full-length near-identical copy (a0); the fragment must NOT inflate it.
        assert!(f1_famcn <= f1_full.len(),
            "famCN bucket (id>=0.98,cov>=0.90) must exclude the half-length fragment: famCN_loci={f1_famcn}, full={} {:?}",
            f1_full.len(), f1);
        assert!(f1.iter().any(|c| c.identity >= 0.98 && c.cov < 0.90),
            "expected the planted fragment to be a >=0.98-identity but <0.90-cov hit (famCN-excluded): {:?}", f1);

        // No cross-family leakage: F1's loci must sit in the A-region of the contig (a0/a8/a15/a_frag
        // offsets, well before the B region begins), and F2's loci must sit in the B region.
        // Layout: 2000 bg | a0 | 20k bg | a8 | 20k bg | a15 | 20k bg | a_frag | 20k bg | b0 | ...
        let b_region_start =
            2_000 + copy_len + 20_000 + copy_len + 20_000 + copy_len + 20_000 + (copy_len / 2) + 20_000;
        assert!(f1.iter().all(|c| c.start < b_region_start),
            "F1 locus leaked into the B region (query-name grouping bug): {:?}", f1);
        assert!(f2.iter().all(|c| c.start >= b_region_start - 100),
            "F2 locus leaked into the A region (query-name grouping bug): {:?}", f2);

        // Bucketing: F1 must have >=2 loci that are >=0.80 identity but <0.98 (the divergent 8%/15%
        // copies) -- proving totalCN (>=0.80) > famCN (>=0.98) is derivable from one projection result.
        let f1_divergent_but_above_floor = f1.iter().filter(|c| c.identity >= 0.80 && c.identity < 0.98).count();
        assert!(f1_divergent_but_above_floor >= 2,
            "F1 expected >=2 loci with 0.80<=identity<0.98 (totalCN>famCN bucketing), got {}: {:?}",
            f1_divergent_but_above_floor, f1);

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// TDD for the `gw_family_catalog --enumerate-copies` totalCN cov floor: a coverage sweep on real
    /// families showed `min_cov=0.50` INFLATES totalCN with partial/domain-fragment hits (GWFAM18: 16 at
    /// cov0.50 vs 12=truth at cov>=0.70), while `min_cov=0.90` is too strict and drops divergent
    /// full-length copies (GSTM 4 vs 19 truth). `min_cov=0.80` is the sweet spot. This test proves that
    /// policy at the `project_families_batch` call level: with `min_cov=0.80`, the planted ~50%-length
    /// near-identical fragment (cov~0.5) must be EXCLUDED from the batch result entirely, while the
    /// full-length divergent copies (0%/~8%/~15% divergence, cov>=0.80) are retained.
    #[test]
    fn project_families_batch_cov80_excludes_fragments() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }

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
        let a0 = gen_seq(1, copy_len);
        let a8 = mutate(&a0, 0.08, 100);
        let a15 = mutate(&a0, 0.15, 200);
        // Same ~50%-length near-identical fragment as the sibling test: clears the old cov>=0.50 floor
        // but must NOT clear the new cov>=0.80 floor.
        let a_frag = mutate(&a0[..(copy_len as usize / 2)], 0.01, 400);

        let dir = std::env::temp_dir().join(format!("rustle_proj_cov80_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut genome: Vec<u8> = Vec::new();
        genome.extend_from_slice(&gen_seq(9001, 2_000));
        genome.extend_from_slice(&a0);
        genome.extend_from_slice(&gen_seq(9002, 20_000));
        genome.extend_from_slice(&a8);
        genome.extend_from_slice(&gen_seq(9003, 20_000));
        genome.extend_from_slice(&a15);
        genome.extend_from_slice(&gen_seq(9007, 20_000));
        genome.extend_from_slice(&a_frag);
        genome.extend_from_slice(&gen_seq(9004, 2_000));

        let fa = dir.join("g.fa");
        std::fs::write(&fa, format!(">chr1\n{}\n", String::from_utf8(genome).unwrap())).unwrap();

        let consensuses: Vec<(String, Vec<u8>)> = vec![("F1".to_string(), a0.clone())];
        let known: std::collections::HashMap<String, Vec<(String, u64, u64)>> = std::collections::HashMap::new();
        let result = project_families_batch(&consensuses, fa.to_str().unwrap(), &known, 0.80, 0.80, "minimap2", 2).unwrap();

        let f1 = result.get("F1").expect("F1 must be present in batch result");

        // The full-length divergent copies (0%/~8%/~15%) all clear cov>=0.80 and must be present.
        assert!(f1.len() >= 3,
            "expected >=3 full-length loci (0%/~8%/~15% divergent) at min_cov=0.80, got {}: {:?}", f1.len(), f1);
        assert!(f1.iter().all(|c| c.cov >= 0.80),
            "min_cov=0.80 call must not return any sub-0.80-coverage hit: {:?}", f1);
        // The ~50%-length fragment must be gone entirely -- no hit with cov in the fragment's ~0.5 band.
        assert!(f1.iter().all(|c| c.cov < 0.40 || c.cov >= 0.75),
            "min_cov=0.80 must exclude the half-length (cov~0.5) fragment hit: {:?}", f1);

        let _ = std::fs::remove_dir_all(&dir);
    }
}
