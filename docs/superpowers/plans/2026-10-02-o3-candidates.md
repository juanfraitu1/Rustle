# O3 candidate copies (`o3_candidates`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put the validated reference-absent-copy chain into `tools/rustle_pipeline.sh` as a `candidates` stage between O1 and O2, built in Rust with no IsoCon.

**Architecture:** A new library module `src/rustle/vg_family/o3_candidates.rs` (pure functions: k-mer sketches, PAF/cs parsing, clustering, consensus, significance merge, flag/link/merge, union) and a thin binary `src/bin/o3_candidates.rs` (BAM passes, batched minimap2 calls through `run_cache`, outputs). `copy_assign` gets a chromosome-aware overlap; the driver gets the stage, the augmentation, a patch realignment and an `assign` that finally consumes `P.fam.copies.*`.

**Tech Stack:** Rust (anyhow, noodles-bam/sam/core, run_cache), minimap2 2.30 via `RUSTLE_MINIMAP2`, samtools in the driver, Python 3 for the fixture generator.

**Spec:** `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`

## Global Constraints

- Build only with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2`, always `--release` for runs, `cargo test` output captured to a file (a pipe to `tail` hides the exit code). Heavy commands (cargo build/test, minimap2, any run > 3 min or > 2 GB) go through `bash tools/rlock.sh heavy ...`, in the foreground; no `pkill -f`, no background waiters.
- Determinism: seeded sampling (seed 1 over sorted read names), stable orders everywhere; threads only inside minimap2.
- The link cut is `--delta`, default `0.00958`; the merge rule is "coverage of the shorter >= 0.5 and `de` <= delta"; the flag floor is `--min-support 6` reads over a component's clusters (the registered >= 2-transcript floor translated: 2 x IsoCon's 3-read transcript minimum; our clusters merge a copy's isoforms, so a cluster count cannot play that role — ruling R1 in the ledger); clusters need `--min-cluster 3` reads; nets are capped at `--max-reads 1000`. These are the spec's values and the prereg's; do not tune them.
- `copy_assign`'s existing fixtures must stay byte-identical (Task 3).
- The canonical checkout `/mnt/c/Users/jfris/Desktop/Rustle` is never edited; work in `/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto` (branch `machine2/soto-evidence`); commits carry the attribution lines.
- Plan rulings against the spec (same observable behaviour, recorded here): alignment engine = batched minimap2 (spec §5.3/5.5/5.7 said poasta): `de`, coverage and mismatch columns come from minimap2's PAF/`cs`, exactly as in Amendments 7-8; the chromosome check in `copy_assign` is a parallel `read_chroms` slice (spec §6 said a field on `AlignedRead`; 41 literal constructors would change for nothing).

## Review Focus

1. A family whose copies lie on two chromosomes with numerically overlapping coordinates: a read on chromosome A must never be attributed to a copy on chromosome B (Task 3 test `cross_chrom_overlap_is_not_an_overlap`).
2. A net above the cap: the sampled set must be the same across runs and machines (Task 5 test `sampling_is_deterministic`).
3. An unmapped read of 69 bp (the real median): never attributed, never counted (Task 4 test `short_unmapped_reads_are_ignored`).
4. A cluster whose consensus links to a reference locus must not become a candidate even if it is in a component with a new-copy cluster: linked clusters are removed before merging (Task 7 test `linked_cluster_never_enters_a_component`).
5. A candidate name colliding with a FASTA sequence name: the augmentation refuses (Task 9, `augment.sh` exits 2 with the name).

---

### Task 1: Prereg Amendment 12 (acceptance), written before any code runs on the held-out

**Files:**
- Modify: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` (append)

- [ ] **Step 1: Append Amendment 12**

```markdown
## Amendment 12 (2026-10-02): `o3_candidates` (the in-house, IsoCon-free stage) on the same held-out (written before the stage ran on it)

Spec `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`. Substrate: Amendment 7's 53 families, masked genome, the same
59,013 scored reads (`linktest/scored.fa`) given to the stage as the BAM `linktest/R.bam` (reads on the masked genome) with a copies
table made from `panel.json`'s surviving copies (clean intervals, `n_reads` from the BAM) in `P.fam.copies.tsv` format, and the
masked splice index. Arm M = masked genome + `P.cand.contigs.fa` (one union per flagged candidate), components as loci, scored by
`merge_test.py score` semantics (D right / wrong / unplaced; S false moves) with the contigs' D/S labels from their best unmasked hit.
- Flag = a component with >= 6 supporting reads over its clusters (`--min-support 6`, the >= 2-transcript floor translated: 2 x IsoCon's
  3-read minimum; the stage's clusters merge a copy's isoforms, so cluster counts cannot play that role). The >= 2-cluster count is reported beside.
- **A12-1 (adopt):** D right >= 80% of IsoCon's 12,787 (>= 10,230) AND false moves <= 5% of S reads.
- **A12-2 (representative):** >= 95% of the reads of each flagged component's clusters keep an AS on the union >= 0.98 x their best
  AS over the component's cluster consensuses (the `rep_choice.py` measure, run on the stage's own alignments).
- **A12-3 (cost):** wall time of the stage on the 53 families <= 40 min (2 x IsoCon's ~20 min) on this machine, 4 threads.
- Reported: candidates per family, clusters per candidate, the deleted copies with no candidate by cause (no reads in the net / clusters
  below the floor / linked to a survivor), the same numbers at delta/2 and 2 x delta, and the flag counts under the alternative >= 2-cluster floor.
```

- [ ] **Step 2: Commit**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && git add docs/PREREG_rna_allele_haplotype_count_2026-10-01.md && git -c user.name="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.name)" -c user.email="$(git -C /mnt/c/Users/jfris/Desktop/Rustle config user.email)" commit -q -m "Prereg Amendment 12: o3_candidates acceptance on the 53-family held-out (rules fixed before the stage exists)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb"
```

---

### Task 2: Module skeleton, k-mer sketches and unmapped-read attribution

**Files:**
- Create: `src/rustle/vg_family/o3_candidates.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod o3_candidates; // O3 candidate copies: read net -> clusters -> consensus -> flag/link/merge -> union (spec 2026-10-02)`)

**Interfaces:**
- Produces: `pub fn canonical_kmers(seq: &[u8], k: usize) -> Vec<u64>`, `pub fn minimizer_sketch(seq: &[u8], k: usize, w: usize) -> Vec<u64>`, `pub fn sketch_share(a: &[u64], b: &[u64]) -> f64` (shared / len of the shorter sketch), `pub struct FamilyKmerIndex`, `FamilyKmerIndex::build(copies: &[(String, Vec<u8>)], k: usize, max_families: usize) -> Self`, `FamilyKmerIndex::attribute(&self, read: &[u8]) -> Option<String>`.

- [ ] **Step 1: Write the failing tests** (bottom of the new file, `#[cfg(test)] mod tests`)

```rust
#[cfg(test)]
mod tests {
    use super::*;
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed; (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x % 4) as usize] }).collect()
    }
    fn mutate(s: &[u8], rate_per_kb: usize, seed: u64) -> Vec<u8> {
        let mut v = s.to_vec(); let mut x = seed.max(1);
        for i in (0..v.len()).step_by(1000 / rate_per_kb.max(1)) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; v[i] = b"ACGT"[((v[i] as u64 + 1 + x % 3) % 4) as usize]; }
        v
    }
    #[test]
    fn canonical_kmers_are_strand_free() {
        let s = rand_seq(500, 7);
        let rc = crate::vg_family::seq_utils::reverse_complement(&s);
        let mut a = canonical_kmers(&s, 31); let mut b = canonical_kmers(&rc, 31);
        a.sort(); b.sort();
        assert_eq!(a, b);
    }
    #[test]
    fn sketch_share_separates_copies_from_errors() {
        let copy_a = rand_seq(3000, 1);
        let same_copy_read = mutate(&copy_a, 2, 3);        // 0.2%: HiFi-like
        let copy_b = mutate(&copy_a, 20, 5);               // 2% diverged paralog
        let sa = minimizer_sketch(&copy_a, 31, 5);
        assert!(sketch_share(&sa, &minimizer_sketch(&same_copy_read, 31, 5)) > 0.8);
        assert!(sketch_share(&sa, &minimizer_sketch(&copy_b, 31, 5)) < 0.6);
    }
    #[test]
    fn attribution_needs_a_clear_winner_and_ignores_short_reads() {
        let fam1 = rand_seq(2000, 11); let fam2 = rand_seq(2000, 12);
        let idx = FamilyKmerIndex::build(&[("F1".into(), fam1.clone()), ("F2".into(), fam2.clone())], 31, 8);
        assert_eq!(idx.attribute(&mutate(&fam1, 2, 2)).as_deref(), Some("F1"));
        assert_eq!(idx.attribute(&rand_seq(1500, 99)), None);               // nothing hits
        let half: Vec<u8> = fam1[..1000].iter().chain(fam2[..1000].iter()).copied().collect();
        assert_eq!(idx.attribute(&half), None);                             // no 2x winner
        assert_eq!(idx.attribute(&fam1[..69]), None);                       // short_unmapped_reads_are_ignored
    }
}
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2 bash tools/rlock.sh heavy cargo test --release --lib o3_candidates > /tmp/claude-1000/-mnt-c-Users-jfris-Desktop/774b64db-68a6-4b08-af48-b829dee19664/scratchpad/t2.out 2>&1; tail -5 .../t2.out`
Expected: compile error (module missing / functions undefined).

- [ ] **Step 3: Implement**

```rust
//! O3 candidate copies (spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md): the reference-absent-copy chain of
//! PREREG_rna_allele_haplotype_count Amendments 7-11 without IsoCon. Pure functions here; BAM passes and minimap2 calls in the binary.
use std::collections::{HashMap, HashSet};

const KMER_K: usize = 31;
const SKETCH_W: usize = 5;
pub const MIN_UNMAPPED_LEN: usize = 300;
const ATTRIB_MIN_FRAC: f64 = 0.30;
const ATTRIB_MIN_RATIO: f64 = 2.0;

fn code(b: u8) -> Option<u64> { match b { b'A' | b'a' => Some(0), b'C' | b'c' => Some(1), b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None } }

/// Canonical (min of forward / reverse-complement) 2-bit k-mers, k <= 32; windows with N are skipped.
pub fn canonical_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k) + 1);
    for &b in seq {
        match code(b) {
            Some(c) => { fw = ((fw << 2) | c) & mask; rv = (rv >> 2) | ((3 - c) << (2 * (k - 1))); valid += 1; }
            None => { valid = 0; fw = 0; rv = 0; }
        }
        if valid >= k { out.push(fw.min(rv)); }
    }
    out
}

fn mix(x: u64) -> u64 { let mut z = x.wrapping_add(0x9E3779B97F4A7C15); z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }

/// (k, w) minimizers over the canonical k-mers: the smallest hashed k-mer of every window of w consecutive k-mers, deduplicated
/// consecutively; sorted and deduplicated so sketches compare by merge.
pub fn minimizer_sketch(seq: &[u8], k: usize, w: usize) -> Vec<u64> {
    let km: Vec<u64> = canonical_kmers(seq, k).into_iter().map(mix).collect();
    let mut out = Vec::new();
    if km.len() < w { out.extend(km.iter().copied().min()); }
    for win in km.windows(w) { let m = *win.iter().min().unwrap(); if out.last() != Some(&m) { out.push(m); } }
    out.sort_unstable(); out.dedup(); out
}

/// Shared minimizers / size of the smaller sketch (both sorted, deduplicated).
pub fn sketch_share(a: &[u64], b: &[u64]) -> f64 {
    let (mut i, mut j, mut shared) = (0, 0, 0usize);
    while i < a.len() && j < b.len() { if a[i] == b[j] { shared += 1; i += 1; j += 1; } else if a[i] < b[j] { i += 1; } else { j += 1; } }
    shared as f64 / a.len().min(b.len()).max(1) as f64
}

/// Canonical 31-mers of every copy sequence -> the families carrying them (k-mers in > max_families families dropped as repeats).
pub struct FamilyKmerIndex { map: HashMap<u64, Vec<String>>, k: usize }
impl FamilyKmerIndex {
    pub fn build(copies: &[(String, Vec<u8>)], k: usize, max_families: usize) -> Self {
        let mut map: HashMap<u64, Vec<String>> = HashMap::new();
        for (fam, seq) in copies {
            for km in canonical_kmers(seq, k) { let v = map.entry(km).or_default(); if !v.contains(fam) { v.push(fam.clone()); } }
        }
        map.retain(|_, v| v.len() <= max_families);
        FamilyKmerIndex { map, k }
    }
    /// The family with the most k-mer hits when >= 30% of the read's k-mers hit it and it leads the runner-up >= 2x; reads < 300 bp: None.
    pub fn attribute(&self, read: &[u8]) -> Option<String> {
        if read.len() < MIN_UNMAPPED_LEN { return None; }
        let kms = canonical_kmers(read, self.k);
        let mut hits: HashMap<&str, usize> = HashMap::new();
        for km in &kms { if let Some(f) = self.map.get(km) { for fam in f { *hits.entry(fam).or_insert(0) += 1; } } }
        let mut v: Vec<(&str, usize)> = hits.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
        let (best, n) = *v.first()?;
        let second = v.get(1).map(|x| x.1).unwrap_or(0);
        if (n as f64) < ATTRIB_MIN_FRAC * kms.len() as f64 { return None; }
        if second > 0 && (n as f64) < ATTRIB_MIN_RATIO * second as f64 { return None; }
        Some(best.to_string())
    }
}
```

- [ ] **Step 4: Run the tests to see them pass** — same command; expected `test result: ok` for the 3 tests.
- [ ] **Step 5: Commit** — `git add src/rustle/vg_family/o3_candidates.rs src/rustle/vg_family/mod.rs`, message `o3_candidates: module skeleton, canonical k-mers, minimizer sketches, family k-mer index (spec 2026-10-02 §5.2-5.3)` + attribution lines.

---

### Task 3: `copy_assign` chromosome-aware overlap

**Files:**
- Modify: `src/rustle/vg_family/copy_assign_pipeline.rs` (`assign_family_detailed_once` :2223, `best_overlap_copy` :1588, `has_block_in_any_copy` :1564, the call sites :2277, :2395, :2495, :2970)
- Modify: the caller(s) of `assign_family_detailed_once` so they pass the reads' chromosomes (they hold `BamRead`s, which carry `chrom`; find them with `grep -n "assign_family_detailed_once(" src`)

**Interfaces:**
- `assign_family_detailed_once(copies, reads, p, genome, mol_names, read_chroms: Option<&[String]>)`; when `Some`, index `i` is read `i`'s chromosome. `best_overlap_copy_on(read, copies, chrom: Option<&str>)` and `has_block_in_any_copy_on(read, copies, chrom: Option<&str>)` restrict `copies` to those with `c.chrom == chrom` when `chrom` is `Some`; the old two-argument functions stay as wrappers passing `None` (byte-identical behaviour for every existing caller).

- [ ] **Step 1: Write the failing test** (in the existing `#[cfg(test)]` module of `copy_assign_pipeline.rs`; build two `DenovoTranscript`s the way neighbouring tests do, one on `"c1"` at 100-500, one on `"c2"` at 120-520, and an `AlignedRead { ref_start: 150, cigar: vec![('M', 200)], seq: vec![b'A'; 200], qual: vec![] }`)

```rust
#[test]
fn cross_chrom_overlap_is_not_an_overlap() {
    let (a, b) = (mk_copy("c1", 100, 500), mk_copy("c2", 120, 520));   // mk_copy = this module's existing transcript helper
    let read = AlignedRead { ref_start: 150, cigar: vec![('M', 200)], seq: vec![b'A'; 200], qual: vec![] };
    assert_eq!(best_overlap_copy_on(&read, &[&a, &b], Some("c2")), Some(1));
    assert_eq!(best_overlap_copy_on(&read, &[&a, &b], Some("c3")), None);
    assert_eq!(best_overlap_copy_on(&read, &[&a, &b], None), best_overlap_copy(&read, &[&a, &b]));
}
```

- [ ] **Step 2: Run** `cargo test --release --lib cross_chrom_overlap` (captured) — expected: `best_overlap_copy_on` not found.
- [ ] **Step 3: Implement**

```rust
pub(crate) fn best_overlap_copy_on(read: &AlignedRead, copies: &[&DenovoTranscript], chrom: Option<&str>) -> Option<usize> {
    let r_end = read_ref_end(read);
    let (mut best, mut best_ov) = (None, 0i64);
    for (ci, c) in copies.iter().enumerate() {
        if chrom.is_some_and(|ch| ch != c.chrom) { continue; }
        let ov = (r_end.min(c.end) as i64) - (read.ref_start.max(c.start) as i64);
        if ov > best_ov { best_ov = ov; best = Some(ci); }
    }
    best
}
pub(crate) fn best_overlap_copy(read: &AlignedRead, copies: &[&DenovoTranscript]) -> Option<usize> { best_overlap_copy_on(read, copies, None) }
```

Same pattern for `has_block_in_any_copy_on`. In `assign_family_detailed_once`, add the parameter `read_chroms: Option<&[String]>` and at each of the four call sites pass `read_chroms.map(|c| c[i].as_str())` (use the read's index there; at :2277 the loop already has it, at :2495 and :2970 it is `read_index` / the enumerate index). Callers: pass `Some(&chroms)` built from the `BamRead`s in the same order as `reads`.

- [ ] **Step 4: Run the whole suite, byte-identical check** — `cargo test --release` captured to a file; expected all green (843 + 1). Then run `tests/copy_assign_families.rs` fixtures and diff their outputs against a pre-change run saved in step 0 (`git stash`-free: run the suite once BEFORE editing and keep the fixture outputs in the scratchpad) — expected: no diff.
- [ ] **Step 5: Commit** — `copy_assign: chromosome-aware overlap (read_chroms slice); byte-identical on single-chromosome families (spec §6)`.

---

### Task 4: PAF and `cs` parsing, pairwise quantities

**Files:**
- Modify: `src/rustle/vg_family/o3_candidates.rs`

**Interfaces:**
- `pub struct PafHit { pub q: String, pub qlen: usize, pub qs: usize, pub qe: usize, pub strand: u8, pub t: String, pub tlen: usize, pub ts: usize, pub te: usize, pub matches: usize, pub block: usize, pub de: f64, pub cs: Option<String> }`
- `pub fn parse_paf(text: &str) -> Vec<PafHit>`; `pub fn best_by_matches(hits: &[PafHit]) -> HashMap<String, PafHit>` (per query, max matches); `pub fn id_cov(h: &PafHit) -> f64` (= matches/block x (qe-qs)/qlen); `pub fn whole_length_d(h: &PafHit) -> f64` (= 1 - matches/qlen); `pub fn shorter_cov(h: &PafHit) -> f64` (span on the shorter of q/t over its length).
- `pub enum CsOp { Eq(usize), Sub(u8, u8), Ins(Vec<u8>), Del(Vec<u8>), Intron(usize) }`; `pub fn parse_cs(cs: &str) -> Vec<CsOp>`.

- [ ] **Step 1: Failing tests**

```rust
#[test]
fn paf_quantities_match_the_chain() {
    let line = "r1\t1000\t10\t990\t+\tc1\t5000\t100\t1100\t960\t1000\t60\tNM:i:40\tde:f:0.0150\tcs:Z::500*ac:300+tt:179-g:0";
    let h = &parse_paf(line)[0];
    assert!((whole_length_d(h) - 0.04).abs() < 1e-9);          // 1 - 960/1000
    assert!((id_cov(h) - 0.96 * 0.98).abs() < 1e-9);
    assert!((h.de - 0.015).abs() < 1e-9);
    let ops = parse_cs(h.cs.as_deref().unwrap());
    assert!(matches!(ops[1], CsOp::Sub(b'a', b'c')));
    assert!(matches!(&ops[3], CsOp::Ins(v) if v == b"tt"));
    assert!(matches!(&ops[5], CsOp::Del(v) if v == b"g"));
}
```

- [ ] **Step 2: Run** (fail: undefined). **Step 3: Implement** — straightforward tab split (12 fixed columns + tags `de:f:`, `cs:Z:`); `parse_cs` walks `:N`, `*xy`, `+seq`, `-seq`, `~xxNNyy`. **Step 4: Run** (pass). **Step 5: Commit** — `o3_candidates: PAF and cs parsing with the chain's quantities (d, identity x coverage, de)`.

---

### Task 5: Clustering, consensus, significance merge (pure functions over PAF/cs)

**Files:**
- Modify: `src/rustle/vg_family/o3_candidates.rs`

**Interfaces:**
- `pub fn cluster_reads(names: &[String], ava: &[PafHit], delta: f64) -> Vec<Vec<usize>>` — union-find over pairs with `shorter_cov >= 0.5 && de <= delta` (best hit per unordered pair by matches); returns clusters as index lists, each sorted, clusters ordered by (size desc, first index).
- `pub fn consensus_from_template(template: &[u8], member_hits: &[(&[u8], &PafHit)]) -> Vec<u8>` — template-and-vote (spec §5.4) from each member's `cs` against the template (`+` strand hits only; members aligned on `-` are reverse-complemented by the caller before alignment, so all hits are `+`).
- `pub fn distinguishing_columns(cs_a_vs_b: &[CsOp]) -> usize` (substitutions only; indels ignored).
- `pub fn variant_is_real(n_small: usize, n_large: usize, k: usize, eps: f64, alpha: f64) -> bool` — p = P(X >= n_small), X ~ Binomial(n_small + n_large, eps^k); real iff p < alpha (k = 0 -> false).
- `pub fn refine_cluster(template: &[u8], members: &[(&[u8], &PafHit)], delta: f64) -> (Vec<usize>, Vec<usize>)` — indices whose `de` to the consensus is <= delta vs the rest (one split pass; the rest become a new cluster if >= 3).

- [ ] **Step 1: Failing tests**

```rust
#[test]
fn clustering_splits_two_percent_and_joins_point_two_percent() {
    // 6 reads: r0..r2 from copy A (0.2% errors), r3..r5 from copy B (2% from A); hits written as the minimap2 ava would give them
    let names: Vec<String> = (0..6).map(|i| format!("r{i}")).collect();
    let hit = |a: usize, b: usize, de: f64| PafHit { q: format!("r{a}"), qlen: 3000, qs: 0, qe: 3000, strand: b'+', t: format!("r{b}"), tlen: 3000, ts: 0, te: 3000, matches: 2900, block: 3000, de, cs: None };
    let ava = vec![hit(0,1,0.002), hit(1,2,0.003), hit(0,2,0.002), hit(3,4,0.002), hit(4,5,0.002), hit(0,3,0.021), hit(2,5,0.019)];
    let cl = cluster_reads(&names, &ava, 0.00958);
    assert_eq!(cl, vec![vec![0,1,2], vec![3,4,5]]);
}
#[test]
fn consensus_vote_fixes_errors_and_keeps_majority_indels() {
    let template = b"ACGTACGTACGTTTTTACGTACGT".to_vec();          // template carries a 1-base error at index 4 (A instead of G)
    let truth    = b"ACGTGCGTACGTTTTTACGTACGT".to_vec();
    // three members equal to truth: cs vs template = ":4*ag:19"  (4 eq, sub a->g, 19 eq)
    let h = PafHit { q: "m".into(), qlen: 24, qs: 0, qe: 24, strand: b'+', t: "t".into(), tlen: 24, ts: 0, te: 24, matches: 23, block: 24, de: 0.04, cs: Some(":4*ag:19".into()) };
    let members: Vec<(&[u8], &PafHit)> = vec![(&truth[..], &h), (&truth[..], &h), (&truth[..], &h)];
    assert_eq!(consensus_from_template(&template, &members), truth);
}
#[test]
fn significance_keeps_supported_variants_and_merges_singletons() {
    assert!(variant_is_real(5, 200, 3, 0.001, 0.01));     // 5 reads carrying 3 distinguishing bases: not error
    assert!(!variant_is_real(2, 200, 1, 0.001, 0.01));    // 2 reads, 1 column: P(X>=2) with p=0.001, n=202 ~ 0.018 >= 0.01 -> merge
    assert!(!variant_is_real(3, 10, 0, 0.001, 0.01));
}
```

- [ ] **Step 2: Run** (fail). **Step 3: Implement**

```rust
pub fn shorter_cov(h: &PafHit) -> f64 { if h.qlen <= h.tlen { (h.qe - h.qs) as f64 / h.qlen.max(1) as f64 } else { (h.te - h.ts) as f64 / h.tlen.max(1) as f64 } }

pub fn cluster_reads(names: &[String], ava: &[PafHit], delta: f64) -> Vec<Vec<usize>> {
    let idx: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut best: HashMap<(usize, usize), &PafHit> = HashMap::new();
    for h in ava {
        let (Some(&a), Some(&b)) = (idx.get(h.q.as_str()), idx.get(h.t.as_str())) else { continue };
        if a == b { continue; }
        let key = (a.min(b), a.max(b));
        if best.get(&key).map_or(true, |o| h.matches > o.matches) { best.insert(key, h); }
    }
    let mut par: Vec<usize> = (0..names.len()).collect();
    fn find(p: &mut Vec<usize>, mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
    for ((a, b), h) in &best { if shorter_cov(h) >= 0.5 && h.de <= delta { let (ra, rb) = (find(&mut par, *a), find(&mut par, *b)); if ra != rb { par[ra.max(rb)] = ra.min(rb); } } }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..names.len() { let r = find(&mut par, i); groups.entry(r).or_default().push(i); }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    for g in &mut out { g.sort_unstable(); }
    out.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));
    out
}

/// Template-and-vote: per template column the majority base over covering members (>= 3 covering, else the template base);
/// insertions < 20 bp after a column present in >= 50% of the members covering it are inserted (the most frequent sequence); deletions
/// < 20 bp in >= 50% remove the column. Indels >= 20 bp are STRUCTURE (isoforms): an insertion >= 20 bp carried by >= 3 members is
/// inserted whatever its share (an exon the template lacks); a deletion >= 20 bp is never applied (the template's exon stays), so the
/// cluster consensus is the exon union of its reads' isoforms with SNV-level majority voting (ruling R2). Members' `cs` strings are relative to the template (`+` strand, ts..te on the template).
pub fn consensus_from_template(template: &[u8], member_hits: &[(&[u8], &PafHit)]) -> Vec<u8> {
    let n = template.len();
    let mut cover = vec![0usize; n]; let mut subs: Vec<HashMap<u8, usize>> = vec![HashMap::new(); n];
    let mut dels = vec![0usize; n]; let mut big_del = vec![false; n]; let mut ins: Vec<HashMap<Vec<u8>, usize>> = vec![HashMap::new(); n + 1];
    for (_, h) in member_hits {
        let Some(cs) = h.cs.as_deref() else { continue };
        let mut t = h.ts;
        for op in parse_cs(cs) {
            match op {
                CsOp::Eq(len) => { for p in t..(t + len).min(n) { cover[p] += 1; } t += len; }
                CsOp::Sub(_, qb) => { if t < n { cover[t] += 1; *subs[t].entry(qb.to_ascii_uppercase()).or_insert(0) += 1; } t += 1; }
                CsOp::Del(seq) => { let big = seq.len() >= 20; for p in t..(t + seq.len()).min(n) { cover[p] += 1; if big { big_del[p] = true; } else { dels[p] += 1; } } t += seq.len(); }
                CsOp::Ins(seq) => { *ins[t.min(n)].entry(seq.to_ascii_uppercase()).or_insert(0) += 1; }
                CsOp::Intron(len) => { t += len; }
            }
        }
    }
    let mut out = Vec::with_capacity(n + 64);
    for p in 0..=n {
        if let Some((seq, cnt)) = ins[p].iter().max_by_key(|(s, c)| (**c, std::cmp::Reverse((*s).clone()))) {
            let covering = if p < n { cover[p] } else { cover[n - 1] };
            if (seq.len() >= 20 && *cnt >= 3) || (seq.len() < 20 && covering >= 3 && 2 * cnt >= covering) { out.extend_from_slice(seq); }
        }
        if p == n { break; }
        if cover[p] >= 3 && 2 * dels[p] >= cover[p] && !big_del[p] { continue; }
        let mut base = template[p].to_ascii_uppercase();
        if cover[p] >= 3 {
            let same = cover[p] - subs[p].values().sum::<usize>() - dels[p];
            if let Some((b, c)) = subs[p].iter().max_by_key(|(b, c)| (**c, std::cmp::Reverse(**b))) { if *c > same { base = *b; } }
        }
        out.push(base);
    }
    out
}

pub fn distinguishing_columns(ops: &[CsOp]) -> usize { ops.iter().filter(|o| matches!(o, CsOp::Sub(..))).count() }

/// P(X >= n_small) for X ~ Binomial(n_small + n_large, eps^k), in log space with a log-factorial table; real iff < alpha.
pub fn variant_is_real(n_small: usize, n_large: usize, k: usize, eps: f64, alpha: f64) -> bool {
    if k == 0 || n_small == 0 { return false; }
    let n = n_small + n_large; let p = eps.powi(k as i32); let (lp, lq) = (p.ln(), (1.0 - p).ln());
    let mut lf = vec![0f64; n + 1]; for i in 1..=n { lf[i] = lf[i - 1] + (i as f64).ln(); }
    let tail: f64 = (n_small..=n).map(|x| (lf[n] - lf[x] - lf[n - x] + x as f64 * lp + (n - x) as f64 * lq).exp()).sum();
    tail < alpha
}
```

`refine_cluster`: members whose hit against the consensus has `de > delta` or `shorter_cov < 0.5` go to the second list.

- [ ] **Step 4: Run** (pass). **Step 5: Commit** — `o3_candidates: read clustering at delta, template-and-vote consensus, real-vs-error merge test`.

---

### Task 6: Flag / link / merge / floor and the exon-union representative

**Files:**
- Modify: `src/rustle/vg_family/o3_candidates.rs`

**Interfaces:**
- `pub struct ClusterSeq { pub family: String, pub id: String, pub n_reads: usize, pub seq: Vec<u8> }`
- `pub enum Fate { InReference, Linked { locus: String, d: f64 }, NewCopy { nearest: String, d: f64 } }`; `pub fn classify(c: &ClusterSeq, best_genome_hit: Option<&PafHit>, delta: f64) -> Fate` — `InReference` iff `id_cov >= 0.999`; `Linked` iff `whole_length_d <= delta`; else `NewCopy`.
- `pub fn components(ids: &[String], ava: &[PafHit], delta: f64) -> Vec<Vec<usize>>` — reuse `cluster_reads` (same rule).
- `pub fn is_flagged(component_clusters: &[&ClusterSeq], min_support: usize) -> bool` — the summed `n_reads` of the component's clusters >= min_support (ruling R1; `--min-support 6`).
- `pub fn union_sequence(members_longest_first: &[Vec<u8>], hits_vs_current: impl FnMut(&[u8], &[u8]) -> Option<PafHit>) -> Vec<u8>` — backbone = first; for each next member the closure aligns it to the CURRENT union (the binary runs minimap2 `-c --cs -x splice:hq`? no: `-x asm20 -c --cs`, member as query, union as target); every `Ins` >= 20 bp is spliced into the union at its target position, an unaligned prefix (`qs >= 20`) is prepended, an unaligned suffix (`qlen - qe >= 20`) appended; positions processed from the end so earlier offsets stay valid.

- [ ] **Step 1: Failing tests**

```rust
#[test]
fn linked_cluster_never_enters_a_component() {
    let c = ClusterSeq { family: "F".into(), id: "F:c0".into(), n_reads: 5, seq: vec![b'A'; 1000] };
    let h = PafHit { q: "F:c0".into(), qlen: 1000, qs: 0, qe: 1000, strand: b'+', t: "chr1".into(), tlen: 1_000_000, ts: 10, te: 1010, matches: 995, block: 1000, de: 0.005, cs: None };
    assert!(matches!(classify(&c, Some(&h), 0.00958), Fate::Linked { .. }));
    let h2 = PafHit { matches: 950, de: 0.05, ..h.clone() };
    assert!(matches!(classify(&c, Some(&h2), 0.00958), Fate::NewCopy { .. }));
    assert!(matches!(classify(&c, None, 0.00958), Fate::NewCopy { .. }));
}
#[test]
fn union_contains_each_exon_once() {
    let e1 = b"ACGTACGTACGTACGTACGTACGT".to_vec(); let e2 = b"TTGACCATGACCATGACCATGACC".to_vec(); let e3 = b"GGCATTGGCATTGGCATTGGCATT".to_vec();
    let iso_a: Vec<u8> = [e1.clone(), e3.clone()].concat();                 // skips e2
    let iso_b: Vec<u8> = [e1.clone(), e2.clone(), e3.clone()].concat();
    // the closure plays minimap2: iso_b vs union(iso_a) has a 24-bp insertion after e1
    let u = union_sequence(&[iso_b.clone(), iso_a.clone()], |_m, _u| None);  // longest first: nothing to add from iso_a
    assert_eq!(u, iso_b);
    let u2 = union_sequence(&[iso_a.clone(), iso_b.clone()], |_m, _u| Some(PafHit { q: "b".into(), qlen: 72, qs: 0, qe: 72, strand: b'+', t: "u".into(), tlen: 48, ts: 0, te: 48, matches: 48, block: 72, de: 0.0, cs: Some(format!(":24+{}:24", String::from_utf8_lossy(&e2).to_lowercase())) }));
    assert_eq!(u2, iso_b);
}
```

- [ ] **Step 2: Run** (fail). **Step 3: Implement** `classify`, `components` (delegates to `cluster_reads`), `union_sequence` (walk `cs` with target position; collect `(pos, seq)` insertions >= 20; apply from the highest position; prefix/suffix from `qs` / `qlen - qe`). **Step 4: Run** (pass). **Step 5: Commit** — `o3_candidates: flag/link/merge fates and the exon-union representative`.

---

### Task 7: minimap2 runner with `run_cache`, and the output writers

**Files:**
- Modify: `src/rustle/vg_family/o3_candidates.rs` (runner + writers)
- Modify: `src/rustle/vg_family/run_cache.rs:320` (`"cand" => &["candidates.tsv", "contigs.fa"]`)

**Interfaces:**
- `pub fn minimap2(args: &[&str], target: &Path, query: &Path, out_paf: &Path, cache: Option<&CacheRoot>) -> Result<()>` — `RUSTLE_MINIMAP2`, `std::process::Command`, stderr to null, `ensure!(status.success())`; with a cache root it keys `kind="paf"` on (args, minimap2 version via `rc::minimap2_version`, FNV hash of target + query bytes) exactly as `mcl_families.rs:2691-2741` and replays by hard link.
- Fixed argument sets (constants): `MM2_AVA = ["-x","asm20","-c","--cs","-X","-N","100","-p","0.1","--secondary=yes"]`, `MM2_MEMBERS = ["-x","asm20","-c","--cs","-N","5","-p","0.5"]` (members vs templates / union), `MM2_GENOME = ["-x","splice:hq","-uf","-c","--eqx","-N","20"]` (consensus vs the splice index).
- `pub struct Candidate { pub family: String, pub id: String, pub clusters: Vec<ClusterSeq>, pub union: Vec<u8>, pub flagged: bool, pub nearest: String, pub d: f64 }`
- `pub fn write_outputs(prefix: &str, cands: &[Candidate], linked: &[(ClusterSeq, String, f64)], nets_for_patch: &[(String, Vec<(String, Vec<u8>)>)]) -> Result<()>` writing `P.cand.candidates.tsv` (`family candidate n_clusters n_reads flagged union_len nearest_locus d n_net n_used`; `flagged` = `n_reads >= --min-support`), `P.cand.clusters.tsv` (`family cluster candidate n_reads consensus_len fate linked_to d`), `P.cand.contigs.fa` (flagged only, `>cand_<family>_<k>`), `P.cand.nets.fa`.

- [ ] **Step 1: Failing test** — `write_outputs` on two synthetic candidates produces files whose headers equal the strings above and whose FASTA holds only the flagged one (`tempdir`).
- [ ] **Step 2–5:** implement, run, commit — `o3_candidates: cached minimap2 runner and output writers; run_cache kind cand`.

---

### Task 8: The binary, the fixture and the integration test

**Files:**
- Create: `src/bin/o3_candidates.rs`
- Modify: `Cargo.toml` (`[[bin]] name = "o3_candidates" path = "src/bin/o3_candidates.rs"`)
- Create: `tests/fixtures/o3_candidates/make_fixture.py`, the generated `genome.fa(.fai)`, `reads.bam(.bai)`, `copies.tsv`, `copies.fa`, `genome.splice.mmi` is NOT committed (built in the test from `genome.fa` with `minimap2 -x splice -d` into a tempdir — the fixture genome is ~60 kb)
- Create: `tests/o3_candidates.rs`

**Binary flow (`main`):** parse args (hand-rolled as `missing_copy_flag.rs:49-66`: `--bam --fasta --copies --copies-fa --index --out [--delta 0.00958] [--max-reads 1000] [--min-cluster 3] [--min-support 6] [--threads 4] [--families]`), load copies (`catalog_input.rs` parsers, column names), per family intervals (`locus_start/locus_end` else `start/end`); BAM pass A (`crate::bam::open_bam`, `reader.query(&header, &index, &region)` per interval, `aligned_read_from_record`; primaries give `seq` oriented as sequenced: reverse-complement when the record's flags say reverse); BAM pass B (one `reader.record_bufs(&header)` sweep: sequences of secondary-only names + unmapped records >= 300 bp through `FamilyKmerIndex::attribute`); cap (sort names, seeded shuffle, `--max-reads`); per family: write `net.fa` -> `minimap2(MM2_AVA, net, net)` -> `cluster_reads` -> drop clusters < `--min-cluster` -> template = longest member -> `minimap2(MM2_MEMBERS, templates.fa, net.fa)` -> `consensus_from_template` per cluster -> `refine_cluster` once -> consensus ava (`MM2_AVA`) -> `variant_is_real` merges (eps 0.001, alpha = the gate's constant; re-polish merged clusters on the larger template) -> all consensuses of all families to one FASTA -> `minimap2(MM2_GENOME, index, consensus.fa)` -> `classify` -> per family `components` over the `NewCopy` clusters (`MM2_AVA` of the family's new-copy consensuses) -> flagged iff the component's clusters sum to >= `--min-support` reads (`is_flagged`) -> union (`MM2_MEMBERS` member vs current union, iteratively) -> `write_outputs`. Temporary files live under `<out>.tmp/` and are removed on success; the whole result is a `run_cache` entry of kind `cand`.

- [ ] **Step 1: Fixture generator** (`make_fixture.py`, run once, outputs committed): a 60 kb random genome with a 3-exon gene (exons 300/200/400 bp) at 10 kb (copy A) and a second copy at 40 kb diverged 3% in the exons (copy B); `genome.fa` keeps only copy A (copy B's region is written as random sequence); 60 reads per copy = spliced transcripts with 0.2% random substitutions and varying 5' starts; reads aligned with `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` to `genome.fa`, sorted, indexed; `copies.tsv` with one family `MCL0`, one copy (A) in the `P.fam.copies.tsv` column layout (`mcl_families.rs:764`), `copies.fa` with A's spliced exon sum. Record the python/minimap2/samtools versions in a `README` line.
- [ ] **Step 2: Failing integration test**

```rust
use std::process::Command;
#[test]
fn flags_the_deleted_copy_and_assign_places_its_reads() {
    let dir = tempfile::tempdir().unwrap();
    let fx = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/o3_candidates");
    let mmi = dir.path().join("genome.splice.mmi");
    assert!(Command::new(std::env::var("RUSTLE_MINIMAP2").unwrap_or("minimap2".into())).args(["-x", "splice", "-d"]).arg(&mmi).arg(format!("{fx}/genome.fa")).status().unwrap().success());
    let out = dir.path().join("t.cand");
    let o = Command::new(env!("CARGO_BIN_EXE_o3_candidates")).args(["--bam", &format!("{fx}/reads.bam"), "--fasta", &format!("{fx}/genome.fa"), "--copies", &format!("{fx}/copies.tsv"), "--copies-fa", &format!("{fx}/copies.fa"), "--index", mmi.to_str().unwrap(), "--out", out.to_str().unwrap(), "--threads", "2"]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let cands = std::fs::read_to_string(format!("{}.candidates.tsv", out.display())).unwrap();
    let flagged: Vec<&str> = cands.lines().skip(1).filter(|l| l.split('\t').nth(4) == Some("1")).collect();
    assert_eq!(flagged.len(), 1, "{cands}");
    let fa = std::fs::read_to_string(format!("{}.contigs.fa", out.display())).unwrap();
    assert!(fa.starts_with(">cand_MCL0_0\n"));
    let union_len = fa.lines().nth(1).unwrap().len();
    assert!((850..=950).contains(&union_len), "union {union_len} bp, expected the 900-bp spliced copy");
}
```

- [ ] **Step 3: Run** (fail: binary missing). **Step 4: Implement the binary** as described. **Step 5: Run** the test (pass) and `cargo test --release` (all green, captured). **Step 6: Commit** — `o3_candidates binary + fixture + integration test (one deleted copy flagged, union ~900 bp)`.

---

### Task 9: Driver — `assign` on `P.fam.copies.*`, the `candidates` stage, augmentation, patch realignment, `flag` corroboration

**Files:**
- Modify: `tools/rustle_pipeline.sh` (flags :70-81, `stage_assign` :272-278, new `stage_candidates`, `stage_flag` :279-291, the `case` :292-296)
- Create: `tools/o3_augment.py` (augmentation: copies rows, copies.fa entries, regions, name-collision check)
- Modify: `src/bin/missing_copy_flag.rs` (`--candidates P.cand.candidates.tsv` -> extra column `o3_candidate`)

- [ ] **Step 1: `stage_assign` reads the default O1 output**

```bash
stage_assign() {
  local copies_tsv="$OUT.fam.copies.tsv" copies_fa="$OUT.fam.copies.fa" regions="$OUT.regions.txt"
  if [ "$LEGACY_CATALOG" = 1 ]; then
    copies_tsv="$OUT.cat.copies.tsv"; copies_fa="$OUT.cat.copies.fa"
    samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":1-"$3}' > "$regions"
  else
    [ -s "$copies_tsv" ] || { echo "assign needs $copies_tsv (run the families stage with a mcl_families that writes the copy table)" >&2; exit 2; }
    cut -f2 "$OUT.fam.copies.regions" > "$regions"      # `{fid}\t{chrom}:{lo}-{hi}`: copy_assign's parse_region takes the FIRST token
  fi
  if [ -s "$OUT.cand.candidates.tsv" ] && [ "$NO_CANDIDATES" != 1 ]; then
    assign_with_candidates "$copies_tsv" "$copies_fa" "$regions"; return
  fi
  say "assign: per-read copy assignment on $copies_tsv"
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$regions" --families "$copies_tsv" --copies-fa "$copies_fa" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l)"
}
```

Flags added to the parser: `--legacy-catalog) LEGACY_CATALOG=1; shift;;`, `--no-candidates) NO_CANDIDATES=1; shift;;`, `--delta) DELTA=$2; shift 2;;` (default `0.00958`), `--cand-max-reads) CAND_MAX=$2; shift 2;;` (default 1000).

- [ ] **Step 2: `stage_candidates` and `assign_with_candidates`**

```bash
stage_candidates() {
  [ "$NO_CANDIDATES" = 1 ] && { say "candidates: skipped (--no-candidates)"; return; }
  [ -n "$INDEX" ] || { echo "candidates needs --index (splice .mmi of the primary genome)" >&2; exit 2; }
  say "candidates: o3_candidates on $OUT.fam.copies.tsv"
  "$BIN/o3_candidates" --bam "$BAM" --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" --index "$INDEX" \
    --delta "$DELTA" --max-reads "$CAND_MAX" --threads "$THREADS" --out "$OUT.cand" > "$OUT.candidates.log" 2>&1
  local n; n=$(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | wc -l)
  say "candidates: $n flagged candidate copies in $(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | cut -f1 | sort -u | wc -l) families"
  [ "$n" -gt 0 ] || return
  python3 "$(dirname "$0")/o3_augment.py" --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" --regions "$OUT.fam.copies.regions" --cand "$OUT.cand" --out "$OUT.aug" || exit $?
  samtools faidx "$OUT.aug.fa"
  say "candidates: realigning the nets of the candidate families to $OUT.aug.fa"
  minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t "$THREADS" "$OUT.aug.fa" "$OUT.cand.nets.fa" 2> "$OUT.aug.mm2.log" | samtools sort -@ 2 -o "$OUT.aug.bam" - && samtools index "$OUT.aug.bam"
}
assign_with_candidates() {
  local copies_tsv=$1 copies_fa=$2 regions=$3
  say "assign: candidate families on $OUT.aug.*, the rest on $copies_tsv"
  "$BIN/copy_assign" --bam "$OUT.aug.bam" --fasta "$OUT.aug.fa" --regions "$OUT.aug.regions.txt" --families "$OUT.aug.copies.tsv" --copies-fa "$OUT.aug.copies.fa" --only-families "$OUT.aug.families.txt" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_cand" > "$OUT.assign_cand.log" 2>&1
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$regions" --families "$copies_tsv" --copies-fa "$copies_fa" --skip-families "$OUT.aug.families.txt" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_rest" > "$OUT.assign_rest.log" 2>&1
  for t in assignments families quant family_join famcn_readonly; do
    { cat "$OUT.assign_cand.$t.tsv"; awk 'NR>1' "$OUT.assign_rest.$t.tsv"; } > "$OUT.assign.$t.tsv"
  done
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows"
}
```

`copy_assign` gains `--only-families FILE` / `--skip-families FILE` (one family id per line) in `load_supplied_families` (`copy_assign.rs:1894`): filter the parsed rows before the contract checks. `o3_augment.py`: reads `P.cand.candidates.tsv` + `P.cand.contigs.fa`; refuses (exit 2, naming it) if a `cand_*` name exists in `FASTA.fai`; writes `P.aug.fa` (= `cat FASTA P.cand.contigs.fa`), `P.aug.copies.tsv` (= copies + one row per flagged candidate: `family_id copy_idx=<max+1..> tid=cand_<f>_<k> chrom=cand_<f>_<k> start=0 end=<len> n_exon=1 strand=+ n_reads=0 exons=0-<len> max_family_identity=0 source=o3_candidate gene_id=. core_hull=. sd_depth=0 core_bp=0 rep_frac=0 member_status=candidate locus_start=0 locus_end=<len>`), `P.aug.copies.fa` (= copies.fa + `>{fid}|{idx}|cand_<f>_<k>:0-<len>|+|nexon=1`), `P.aug.regions.txt` (= second column of `P.fam.copies.regions` + `cand_<f>_<k>:0-<len>`), `P.aug.families.txt`. Verify `member_status=candidate` passes `catalog_input.rs:186/298` (it must be accepted like `partner`: skip the read check; add it there if not).

- [ ] **Step 3: `stage_flag` corroboration** — `missing_copy_flag --candidates P.cand.candidates.tsv`: after the verdict table is built, add column `o3_candidate` = the id of a flagged candidate whose `nearest_locus` equals the row's locus name (else `-`). Unit test in `missing_copy_flag.rs`'s test module on two synthetic rows.
- [ ] **Step 4: Stage wiring** — `case` gets `candidates) stage_candidates;;` and `all) stage_assemble; stage_families; stage_candidates; stage_assign; stage_flag;;`; header comment updated (O1 families -> O3 candidates -> O2 assign -> O3 flag). `bash -n tools/rustle_pipeline.sh`.
- [ ] **Step 5: End-to-end on the fixture** — run `tools/rustle_pipeline.sh --bam tests/fixtures/o3_candidates/reads.bam --fasta tests/fixtures/o3_candidates/genome.fa --index <tmp>/genome.splice.mmi --out <tmp>/fx all --bin $CARGO_TARGET_DIR/release` after `cargo build --release`; expected: `candidates: 1 flagged`, `assign:` rows for `MCL0` include reads assigned to `cand_MCL0_0`. Save the log in the scratchpad.
- [ ] **Step 6: Commit** — `pipeline: assign consumes P.fam.copies.*; candidates stage with augmentation and patch realignment; flag corroboration column; copy_assign --only/--skip-families`.

---

### Task 10: Acceptance (Amendment 12) on the 53-family held-out

**Files:**
- Create: `bench/rna_allele/accept_o3_candidates.sh` (the exact commands below), `bench/rna_allele/panel_to_copies.py`
- Create: `docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`; modify `docs/NEGATIVE_RESULTS_REGISTER.md` (rows 1216+)

- [ ] **Step 1: Inputs** — `panel_to_copies.py` turns `linktest/panel.json`'s surviving copies into `A12.copies.tsv` (`P.fam.copies.tsv` columns; `start/end` = clean interval, `exons` = `start-end`, `n_reads` counted from `linktest/R.bam`, `source=panel`) and `A12.copies.fa` (genomic interval sequence from `linktest/masked.fa`, `+`, `nexon=1`); `A12.regions` as `{fid}\t{chrom}:{lo-5000}-{hi+5000}`.
- [ ] **Step 2: Run the stage** (heavy, foreground, timed): `o3_candidates --bam linktest/R.bam --fasta linktest/masked.fa --copies A12.copies.tsv --copies-fa A12.copies.fa --index linktest/masked.splice.mmi --out A12.cand --threads 4` under `/usr/bin/time -v`; record wall time (A12-3).
- [ ] **Step 3: Arm M** — `cat linktest/masked.fa A12.cand.contigs.fa > A12.M.fa`, splice index, realign `linktest/scored.part{0,1,2}.fa` (three heavy calls, the usual flags), merge to `A12.M.bam`; label every contig (`cand_*`) by its best hit in the UNMASKED genome (`minimap2 -c -x splice:hq -uf -N 20 GGO.splice.mmi A12.cand.contigs.fa`): D-derived iff on the family's masked interval; build a `contigs.tsv` in `link_test.py`'s layout (one row per candidate, `linked=0`) and a one-row-per-candidate `merge/components.tsv`, then `merge_test.py score --w A12dir` (R.bam symlinked) — D right / wrong / unplaced, S false moves.
- [ ] **Step 4: A12-2** — align the clusters' reads (from `A12.cand.clusters.tsv` membership) to the unions and to their cluster consensuses (`minimap2 -c -x asm20`), compute the kept fraction as `rep_choice.py` does.
- [ ] **Step 5: Verdicts and write-up** — A12-1/2/3 as registered; the cause table for deleted copies without a candidate; delta/2 and 2 x delta reruns of the stage (`--delta`) for the reported-beside line; doc, register rows, memory file update, commit, push (branch and `main`).

---

## Self-review notes

- Spec coverage: §5.1-5.2 (Task 2 + binary), §5.3-5.5 (Task 5), §5.6-5.7 (Task 6), §5.8 (Task 7), §6 (Task 3), §7 (Task 9), §8 (binary + augment refusal), §9 (Tasks 2-8 unit/fixture; Task 10 acceptance; Task 1 prereg). §10 deferred items untouched.
- Type consistency: `PafHit` fields used identically in Tasks 4-7; `ClusterSeq`/`Candidate` defined in Task 6/7 and consumed by the binary (Task 8) and `write_outputs` (Task 7).
- The two rulings (minimap2 engine; `read_chroms` slice) are stated in Global Constraints and do not change any registered rule.
