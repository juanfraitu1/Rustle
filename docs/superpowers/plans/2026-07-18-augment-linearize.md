# Augment-and-Linearize Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Emit a threshold-free "linearize certificate" for each candidate reference-absent copy — augment the local family reference with the candidate, re-align the ambiguous MAPQ-0 read pool, and score how uniquely they land on the new copy vs a dinucleotide-shuffled decoy.

**Architecture:** A new pure module `linearize.rs` (deterministic decoy shuffle + a `linearize_certificate` core that takes an injected realign closure) keeps the statistics minimap2-free and unit-testable; `absent_copy.rs` supplies the real minimap2 realign closure; `denovo_pipeline.rs` wires it at the stage-2 admission call site; `copy_assign.rs` writes `<out>.linearize.tsv` and adds the opt-in `--linearize-gate`.

**Tech Stack:** Rust (crate `rustle`), minimap2 (PAF mode), the existing `remap_identity_minimap2`/`parse_paf_hits` shell + temp-FASTA patterns.

## Global Constraints

- Default `copy_assign` (no `--absent-copies`) byte-identical; the certificate runs only in the absent-copy path.
- The certificate core is a PURE function with an injected realign closure (minimap2-free unit tests).
- "Unique/linearized" = a read whose PRIMARY alignment lands on the candidate contig with `MAPQ > 0`.
- Deterministic decoy shuffle (reproducible certificates); seed = a hash of the candidate sequence.
- Report-first by default (`<out>.linearize.tsv`); admission decisions unchanged unless `--linearize-gate`.
- No new k-mer computation; reuse the existing minimap2 shell + PAF-parse infra.

## File Structure

- **Create** `src/rustle/vg_family/linearize.rs` — `Lcg`/`fisher_yates`, `dinucleotide_shuffle`, `LinearizeCertificate`, `Verdict`, `linearize_certificate` (pure core). Add `pub mod linearize;` to `src/rustle/vg_family/mod.rs`.
- **Modify** `src/rustle/vg_family/absent_copy.rs` — `realign_pool_minimap2` (the real closure) + `linearize_certificate_minimap2` wrapper.
- **Modify** `src/rustle/vg_family/denovo_pipeline.rs` — build the pool + call the certificate at the stage-2 call site; thread a `Vec<LinearizeCertificate>` out of `detect_and_assign`.
- **Modify** `src/bin/copy_assign.rs` — `<out>.linearize.tsv` writer; `--linearize-gate` flag.

---

### Task 1: Deterministic decoy generator (dinucleotide shuffle + RC)

**Files:**
- Create: `src/rustle/vg_family/linearize.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod linearize;`)
- Test: in `linearize.rs`

**Interfaces:**
- Produces: `pub fn dinucleotide_shuffle(seq: &[u8], seed: u64) -> Vec<u8>` — same length, same first/last base, same exact dinucleotide counts, scrambled order, deterministic in `seed`. Reuse `crate::vg_family::seq_utils::reverse_complement` for the RC decoy (do not re-implement).

- [ ] **Step 1: Write the failing test**
```rust
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib linearize::tests::dinucleotide_shuffle_preserves 2>&1 | tail -20`
Expected: FAIL — module/function not found.

- [ ] **Step 3: Write minimal implementation**

At the top of `linearize.rs`:
```rust
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
```
Add `pub mod linearize;` to `src/rustle/vg_family/mod.rs`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib linearize::tests::dinucleotide_shuffle_preserves 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/linearize.rs src/rustle/vg_family/mod.rs
git commit -m "feat(linearize): deterministic dinucleotide-preserving decoy shuffle"
```

---

### Task 2: `LinearizeCertificate` + the pure `linearize_certificate` core

**Files:**
- Modify: `src/rustle/vg_family/linearize.rs`
- Test: same file

**Interfaces:**
- Produces:
```rust
pub enum Verdict { Linearizes, Not, Undetermined }
pub struct LinearizeCertificate {
    pub n_pool: usize,
    pub linearized_frac_real: f64,
    pub mean_frac_decoy: f64,
    pub delta: f64,
    pub perm_p: f64,
    pub verdict: Verdict,
}
/// `realign(ref_contigs, reads)` returns, per read, `Some((best_contig_idx, mapq))` or `None` (unmapped).
/// The candidate is always the LAST contig of each ref set built here (index == family_copy_seqs.len()).
pub fn linearize_certificate(
    candidate_seq: &[u8],
    family_copy_seqs: &[Vec<u8>],
    pool_reads: &[Vec<u8>],
    n_decoys: usize,
    seed: u64,
    min_pool: usize,
    alpha: f64,
    realign: impl Fn(&[Vec<u8>], &[Vec<u8>]) -> Vec<Option<(usize, u32)>>,
) -> LinearizeCertificate
```

- [ ] **Step 1: Write the failing test**
```rust
    // a fake realign: a read "belongs" to the candidate iff its bytes equal the candidate contig's bytes;
    // then it maps uniquely (mapq 60) to the candidate index; otherwise it maps to copy 0 with mapq 0 (tied).
    fn fake_realign(refs: &[Vec<u8>], reads: &[Vec<u8>]) -> Vec<Option<(usize, u32)>> {
        let cand_idx = refs.len() - 1;
        reads.iter().map(|r| {
            if r == &refs[cand_idx] { Some((cand_idx, 60)) } else { Some((0, 0)) }
        }).collect()
    }
    #[test]
    fn real_copy_linearizes_decoy_does_not() {
        let cand = b"ACGTACGTTTGGCCAAACGTACGT".to_vec();
        let copies = vec![b"TTTTGGGGCCCCAAAA".to_vec()];
        // 8 pool reads that ARE the candidate (they linearize on it), 2 that are noise
        let mut pool: Vec<Vec<u8>> = (0..8).map(|_| cand.clone()).collect();
        pool.push(b"NNNNNN".to_vec()); pool.push(b"NNNNNN".to_vec());
        let cert = linearize_certificate(&cand, &copies, &pool, 19, 7, 5, 0.05, fake_realign);
        assert!((cert.linearized_frac_real - 0.8).abs() < 1e-9, "8/10 land on candidate");
        assert_eq!(cert.mean_frac_decoy, 0.0, "decoys != candidate bytes -> no read lands on them");
        assert!(cert.delta > 0.5);
        assert!(cert.perm_p <= 1.0 / 20.0 + 1e-9, "no decoy beats real -> perm_p = 1/(N+1)");
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
        let cert = linearize_certificate(&cand, &[], &[cand.clone(); 2], 19, 7, 5, 0.05, fake_realign);
        assert!(matches!(cert.verdict, Verdict::Undetermined));
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib linearize::tests::real_copy_linearizes 2>&1 | tail -20`
Expected: FAIL — types/fn not found.

- [ ] **Step 3: Write minimal implementation**
```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Verdict { Linearizes, Not, Undetermined }
#[derive(Clone, Debug)]
pub struct LinearizeCertificate {
    pub n_pool: usize,
    pub linearized_frac_real: f64,
    pub mean_frac_decoy: f64,
    pub delta: f64,
    pub perm_p: f64,
    pub verdict: Verdict,
}

fn frac_on_candidate(hits: &[Option<(usize, u32)>], cand_idx: usize) -> f64 {
    if hits.is_empty() { return 0.0; }
    let k = hits.iter().filter(|h| matches!(h, Some((i, q)) if *i == cand_idx && *q > 0)).count();
    k as f64 / hits.len() as f64
}

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
        return LinearizeCertificate { n_pool, linearized_frac_real: f64::NAN, mean_frac_decoy: f64::NAN,
            delta: f64::NAN, perm_p: f64::NAN, verdict: Verdict::Undetermined };
    }
    let cand_idx = family_copy_seqs.len(); // candidate is the LAST contig
    let build = |extra: &[u8]| -> Vec<Vec<u8>> {
        let mut v = family_copy_seqs.to_vec(); v.push(extra.to_vec()); v
    };
    let real = frac_on_candidate(&realign(&build(candidate_seq), pool_reads), cand_idx);
    // decoys: N dinucleotide shuffles (distinct seeds) + the reverse-complement as one extra decoy
    let mut decoy_fracs: Vec<f64> = Vec::with_capacity(n_decoys + 1);
    for d in 0..n_decoys {
        let decoy = dinucleotide_shuffle(candidate_seq, seed.wrapping_add(d as u64 + 1));
        decoy_fracs.push(frac_on_candidate(&realign(&build(&decoy), pool_reads), cand_idx));
    }
    let rc = crate::vg_family::seq_utils::reverse_complement(candidate_seq);
    decoy_fracs.push(frac_on_candidate(&realign(&build(&rc), pool_reads), cand_idx));
    let nd = decoy_fracs.len();
    let mean_decoy = decoy_fracs.iter().sum::<f64>() / nd as f64;
    let n_ge = decoy_fracs.iter().filter(|&&d| d >= real).count();
    let perm_p = (n_ge as f64 + 1.0) / (nd as f64 + 1.0);
    let verdict = if perm_p < alpha { Verdict::Linearizes } else { Verdict::Not };
    LinearizeCertificate { n_pool, linearized_frac_real: real, mean_frac_decoy: mean_decoy,
        delta: real - mean_decoy, perm_p, verdict }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib linearize:: 2>&1 | tail -20`
Expected: PASS (all three).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/linearize.rs
git commit -m "feat(linearize): LinearizeCertificate + pure certificate core (real vs decoy perm-p)"
```

---

### Task 3: The real minimap2 realign closure

**Files:**
- Modify: `src/rustle/vg_family/absent_copy.rs`
- Test: `absent_copy.rs` (a PAF-parse unit test; the subprocess path is exercised by the Task 6 data-gated render)

**Interfaces:**
- Produces: `pub fn realign_pool_minimap2(ref_contigs: &[Vec<u8>], reads: &[Vec<u8>]) -> Vec<Option<(usize, u32)>>` — matches the `realign` closure shape from Task 2; writes `ref_contigs` as a temp target FASTA (`>{i}\n{seq}`), `reads` as a temp query FASTA (`>{j}\n{seq}`), runs minimap2, and for each read `j` returns its PRIMARY hit's `(target_contig_index, mapq)` (or `None` if unmapped).

- [ ] **Step 1: Write the failing test**

A pure parse test (no minimap2) of the PAF→per-read-primary reducer that `realign_pool_minimap2` uses:
```rust
#[test]
fn primary_hit_reducer_picks_tp_P_and_reads_mapq() {
    // read 0: two hits, primary (tp:A:P) on contig 2 mapq 60; read 1: single hit contig 0 mapq 0
    let paf = "0\t100\t0\t100\t+\t2\t120\t0\t100\t95\t100\t60\ttp:A:P\n\
               0\t100\t0\t100\t+\t1\t120\t0\t100\t80\t100\t0\ttp:A:S\n\
               1\t100\t0\t100\t+\t0\t120\t0\t100\t90\t100\t0\ttp:A:P\n";
    let hits = super::paf_primary_hits(paf, 2);   // 2 reads
    assert_eq!(hits[0], Some((2, 60)));
    assert_eq!(hits[1], Some((0, 0)));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib absent_copy::tests::primary_hit_reducer 2>&1 | tail -20`
Expected: FAIL — `paf_primary_hits` not found.

- [ ] **Step 3: Write minimal implementation**

Add to `absent_copy.rs` (contig/read headers are the integer index as text):
```rust
/// PAF -> per-read primary hit. `n_reads` reads were written with headers `>0..>n_reads`; targets with `>0..`.
/// Primary = the `tp:A:P` hit (fall back to the largest alignment-block length, PAF col 10). Returns per read
/// `Some((target_contig_index, mapq))` or `None`. PAF is 0-indexed here: f[0]=qname, f[5]=tname, f[11]=mapq.
pub(crate) fn paf_primary_hits(paf: &str, n_reads: usize) -> Vec<Option<(usize, u32)>> {
    // per read: (is_primary, block_len, contig_idx, mapq)
    let mut best: Vec<Option<(bool, u64, usize, u32)>> = vec![None; n_reads];
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let q: usize = match f[0].parse() { Ok(v) => v, Err(_) => continue };
        if q >= n_reads { continue; }
        let t: usize = match f[5].parse() { Ok(v) => v, Err(_) => continue };
        let block: u64 = f[10].parse().unwrap_or(0);
        let mapq: u32 = f[11].parse().unwrap_or(0);
        let is_p = f[12..].iter().any(|x| *x == "tp:A:P");
        let cand = (is_p, block, t, mapq);
        // prefer primary; then larger block
        let take = match best[q] { None => true, Some((p0, b0, _, _)) => (is_p, block) > (p0, b0) };
        if take { best[q] = Some(cand); }
    }
    best.into_iter().map(|o| o.map(|(_, _, t, mq)| (t, mq))).collect()
}

/// Real realign closure: temp target FASTA (ref_contigs) + temp query FASTA (reads) -> minimap2 PAF ->
/// per-read primary (contig_idx, mapq). Mirrors `remap_identity_minimap2`'s temp-file + spawn pattern.
pub fn realign_pool_minimap2(ref_contigs: &[Vec<u8>], reads: &[Vec<u8>]) -> Vec<Option<(usize, u32)>> {
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NONCE: AtomicUsize = AtomicUsize::new(0);
    if ref_contigs.is_empty() || reads.is_empty() { return vec![None; reads.len()]; }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let n = NONCE.fetch_add(1, Ordering::Relaxed);
    let refp = dir.join(format!("rustle_lin_ref_{pid}_{n}.fa"));
    let qp = dir.join(format!("rustle_lin_q_{pid}_{n}.fa"));
    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup { fn drop(&mut self) { let _ = std::fs::remove_file(&self.0); } }
    let (_c1, _c2) = (Cleanup(refp.clone()), Cleanup(qp.clone()));
    let write_fa = |path: &std::path::Path, seqs: &[Vec<u8>]| -> Option<()> {
        let mut fh = std::fs::File::create(path).ok()?;
        for (i, s) in seqs.iter().enumerate() { writeln!(fh, ">{}", i).ok()?; fh.write_all(s).ok()?; fh.write_all(b"\n").ok()?; }
        Some(())
    };
    if write_fa(&refp, ref_contigs).is_none() || write_fa(&qp, reads).is_none() { return vec![None; reads.len()]; }
    // intronless consensus contigs -> map-hifi-style, NOT -x splice. --secondary=no + default -N so a
    // non-unique read gets mapq 0 (minimap2's coin-flip signal). Query is the reads, target is the ref set.
    let out = match std::process::Command::new(&mm2)
        .args(["-c", "-x", "map-hifi", "--secondary=no", "-t", "1"]).arg(&refp).arg(&qp).output() {
        Ok(o) if o.status.success() => o,
        _ => return vec![None; reads.len()],
    };
    paf_primary_hits(&String::from_utf8_lossy(&out.stdout), reads.len())
}
```
(If `map-hifi` is unavailable in the local minimap2, `asm20` is an acceptable substitute — both are intronless presets; the render in Task 6 will confirm.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib absent_copy::tests::primary_hit_reducer 2>&1 | tail -20` (PASS) + `cargo build -p rustle 2>&1 | tail -5` (Finished).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/absent_copy.rs
git commit -m "feat(absent_copy): minimap2 realign-pool closure + PAF primary-hit reducer (mapq)"
```

---

### Task 4: Wire the certificate at the stage-2 admission call site

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (stage-2 block ~1524-1558; `detect_and_assign` return)
- Test: `denovo_pipeline.rs` (a wiring test using the injected fake realign, no minimap2)

**Interfaces:**
- Consumes: `linearize::linearize_certificate`, `absent_copy::realign_pool_minimap2`, the in-scope `region`/`region_mapq`/`all_copies`.
- Produces: `detect_and_assign` also returns `Vec<(String /*family_id*/, LinearizeCertificate, (String,u64,u64) /*candidate locus*/)>`. Compute for each `Admission::Copy(t, _)` on the stage-2 arm; the pool = `region[i].seq` where `region_mapq[i] == 0`; the family copy seqs = each already-admitted copy's spliced seq.

- [ ] **Step 1: Write the failing test**

Because the stage-2 path needs a full family fixture, test the extracted helper instead: add and test `fn family_linearize_cert(t_seq, copy_seqs, pool, realign) -> LinearizeCertificate` that packages the pool build + call (defaults `n_decoys=19, min_pool=5, alpha=0.05, seed = fnv hash of t_seq`):
```rust
#[test]
fn family_linearize_cert_uses_mapq0_pool() {
    use crate::vg_family::linearize::Verdict;
    let cand = b"ACGTACGTTTGGCCAAACGTACGT".to_vec();
    let pool: Vec<Vec<u8>> = (0..8).map(|_| cand.clone()).collect();
    let realign = |refs: &[Vec<u8>], reads: &[Vec<u8>]| {
        let ci = refs.len() - 1;
        reads.iter().map(|r| if r == &refs[ci] { Some((ci, 60u32)) } else { Some((0, 0)) }).collect::<Vec<_>>()
    };
    let cert = family_linearize_cert(&cand, &[b"TTTTGGGGCCCC".to_vec()], &pool, realign);
    assert!(matches!(cert.verdict, Verdict::Linearizes));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib family_linearize_cert_uses 2>&1 | tail -20`
Expected: FAIL — `family_linearize_cert` not found.

- [ ] **Step 3: Write minimal implementation**

Add `family_linearize_cert` (the seed = a simple FNV-1a of `t_seq` for reproducibility), then wire it in the stage-2 block. In the admission loop (`denovo_pipeline.rs:~1527-1540`), on the `Admission::Copy(t, id)` arm, before `admitted.push(t)`:
```rust
    let pool: Vec<Vec<u8>> = region.iter().zip(region_mapq.iter())
        .filter(|(_, &q)| q == 0).map(|(r, _)| r.seq.clone()).collect();
    let copy_seqs: Vec<Vec<u8>> = all_copies.iter().map(|c| c.seq.clone()).collect();
    let cert = family_linearize_cert(&t.seq, &copy_seqs, &pool, absent_copy::realign_pool_minimap2);
    linearize_certs.push((cf.family_id.clone(), cert.clone(), (t.chrom.clone(), t.start, t.end)));
    // (Task 5 adds: if linearize_gate && cert.verdict != Linearizes { dna_needs.push(...); continue; })
```
Declare `let mut linearize_certs: Vec<(String, LinearizeCertificate, (String,u64,u64))> = Vec::new();` before the family loop, and add it to `detect_and_assign`'s return tuple (update the signature + the `copy_assign.rs` call-site destructure; `cargo build` flags the arity).

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib family_linearize_cert_uses 2>&1 | tail -20` (PASS) + `cargo build -p rustle 2>&1 | tail -5` (Finished) + `cargo test -p rustle --lib denovo_pipeline:: 2>&1 | grep "test result"` (no new regressions; the known `family_merge_default…` fixture failure is unrelated).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo): compute linearize certificate per admitted absent copy (report-first)"
```

---

### Task 5: `<out>.linearize.tsv` writer + `--linearize-gate`

**Files:**
- Modify: `src/bin/copy_assign.rs` (destructure the new return; TSV writer; `--linearize-gate` flag)
- Test: `copy_assign.rs` (a TSV-row formatting unit test)

**Interfaces:**
- Consumes: the `Vec<(String, LinearizeCertificate, (String,u64,u64))>` from `detect_and_assign`.
- Produces: `<out>.linearize.tsv` with header `family_id\tchrom\tstart\tend\tn_pool\tlinearized_frac_real\tmean_frac_decoy\tdelta\tperm_p\tverdict`; a `--linearize-gate` bool arg (default false).

- [ ] **Step 1: Write the failing test**
```rust
#[test]
fn linearize_tsv_row_formats() {
    use rustle::vg_family::linearize::{LinearizeCertificate, Verdict};
    let c = LinearizeCertificate { n_pool: 40, linearized_frac_real: 0.82, mean_frac_decoy: 0.01,
        delta: 0.81, perm_p: 0.05, verdict: Verdict::Linearizes };
    let row = linearize_tsv_row("GWFAM1", ("chr9", 100, 200), &c);
    assert_eq!(row, "GWFAM1\tchr9\t100\t200\t40\t0.820\t0.010\t0.810\t0.0500\tLINEARIZES");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --bin copy_assign linearize_tsv_row_formats 2>&1 | tail -20`
Expected: FAIL — `linearize_tsv_row` not found.

- [ ] **Step 3: Write minimal implementation**
```rust
fn verdict_str(v: rustle::vg_family::linearize::Verdict) -> &'static str {
    use rustle::vg_family::linearize::Verdict::*;
    match v { Linearizes => "LINEARIZES", Not => "NOT", Undetermined => "UNDETERMINED" }
}
fn linearize_tsv_row(fam: &str, loc: (&str, u64, u64), c: &rustle::vg_family::linearize::LinearizeCertificate) -> String {
    let f = |x: f64| if x.is_nan() { "NA".to_string() } else { format!("{:.3}", x) };
    format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", fam, loc.0, loc.1, loc.2, c.n_pool,
        f(c.linearized_frac_real), f(c.mean_frac_decoy), f(c.delta),
        if c.perm_p.is_nan() { "NA".to_string() } else { format!("{:.4}", c.perm_p) }, verdict_str(c.verdict))
}
```
Add the `--linearize-gate` arg to `struct Args` (`#[arg(long)] linearize_gate: bool`). After the sweep, write `<out>.linearize.tsv` (header + one `linearize_tsv_row` per cert). Pass `args.linearize_gate` into `detect_and_assign` (via the config it already threads for `absent_copies`) so Task 4's TODO fires: when `linearize_gate` and `cert.verdict != Linearizes`, the candidate becomes a `DnaNeeds` record (`reason = "did not linearize (perm_p >= alpha)"`) instead of an admitted copy.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --bin copy_assign linearize_tsv_row_formats 2>&1 | tail -20` (PASS) + `cargo build --release -p rustle --bin copy_assign 2>&1 | tail -5` (Finished — binary for Task 6).

- [ ] **Step 5: Commit**
```bash
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): <out>.linearize.tsv + opt-in --linearize-gate"
```

---

### Task 6: DSFAM26 data-gated validation + final review

**Files:** manual (data-gated)

- [ ] **Step 1: Run on the DSFAM26 collapsed MHC region with --absent-copies**

FOREGROUND, output under winloci_scratch (WSL2 crash rule). DSFAM26's collapsed copy is the real reference-absent case:
```bash
RUSTLE_SKIP_POA_DIAGNOSTIC=1 ./target/release/copy_assign \
  --bam /home/juanfra/winloci_scratch/GGO_mm.bam --fasta /home/juanfra/winloci_scratch/GGO.fasta \
  --region NC_073229.2:47357766-49074923 --out /home/juanfra/winloci_scratch/dsfam26_lin \
  --absent-copies --max-poa-len 80000
column -t /home/juanfra/winloci_scratch/dsfam26_lin.linearize.tsv
```
Expected: at least one candidate with a `LINEARIZES` verdict (`linearized_frac_real` ≫ `mean_frac_decoy`, small `perm_p`). If no families/candidates form on that region (v1 found DSFAM26 is a gw-catalog family that doesn't reconstruct from a plain `--region`), fall back to any region that admits an absent copy under `--absent-copies` and confirm the certificate is emitted and self-consistent (`perm_p` in `[1/(N+1), 1]`, `delta = real − mean_decoy`). Record what ran in the report — the certificate math is already locked by the Task 2/3 hermetic tests; this step confirms the live minimap2 arm produces a sane certificate on real reads.

- [ ] **Step 2: Final whole-branch review** (dispatch the final reviewer), then finish the branch.

---

## Self-Review

**Spec coverage:**
- Deterministic dinucleotide-shuffle decoy (+RC) → Task 1. ✓
- `LinearizeCertificate` + pure core (real vs decoy, perm_p, verdict, UNDETERMINED on small pool) → Task 2. ✓
- Real minimap2 realign closure (augmented FASTA, map-hifi/asm20, PAF col5/col11/tp:A:P) → Task 3. ✓
- linearized = primary on candidate contig with MAPQ>0 → `frac_on_candidate` (Task 2) + `paf_primary_hits` (Task 3). ✓
- Wire at stage-2 call site; pool from region MAPQ-0; report-first Vec → Task 4. ✓
- `<out>.linearize.tsv` + opt-in `--linearize-gate` (non-LINEARIZES → DnaNeeds) → Task 5. ✓
- Hermetic (injected closure) minimap2-free tests → Tasks 2 & 4. ✓
- Default byte-identical (all in the `--absent-copies` path) → Task 4/5 confinement. ✓
- DSFAM26 data-gated validation → Task 6. ✓

**Placeholder scan:** no TBD/TODO; every code step has complete code; the one deferred item (Task 4's gate TODO) is explicitly completed in Task 5.

**Type consistency:** `LinearizeCertificate`/`Verdict`/`linearize_certificate`/`frac_on_candidate`/`dinucleotide_shuffle`/`realign_pool_minimap2`/`paf_primary_hits`/`family_linearize_cert`/`linearize_tsv_row` are stable across tasks. The `realign` closure shape `Fn(&[Vec<u8>], &[Vec<u8>]) -> Vec<Option<(usize,u32)>>` matches between the core (Task 2), the real closure (Task 3), and the wiring (Task 4). One integration note (Task 4/5): `detect_and_assign`'s return tuple grows by one field — every call site must destructure it (`cargo build` enforces).
