# Assembly-based parCN (`parcn`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship an optional standalone `parcn` binary that projects catalog copy consensuses onto the two phased haplotype assemblies and reports per-copy paralog-specific copy number (parCN), split mat/pat, using deterministic SUN witnesses.

**Architecture:** New `src/rustle/vg_family/parcn.rs` (pure core: private-position/SUN computation, cs-tag base reading, hybrid assignment, tabulation) + a CIGAR/cs-retaining projection variant in `genome_projection.rs` + a thin `src/bin/parcn.rs` orchestrator. Consumes only `copies.fa` + two splice `.mmi`; writes only `<out>.parcn*.tsv`; never touches the RNA-exclusive core.

**Tech Stack:** Rust; `minimap2` (external, splice mode, `--cs`); the in-repo banded Gotoh aligner `banded_msa_pair`; `clap` for the binary.

## Global Constraints

- **RNA-exclusive core untouched:** `parcn` adds NO call site to and changes NO output of `gw_family_catalog` / `copy_assign` / the RNA pipeline. Removing `parcn` leaves the pipeline byte-identical. (Exposing an existing private fn as `pub(crate)` is allowed; changing its behavior is not.)
- **Inputs:** `--copies-fa` (header `>{family_id}|{copy_idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`), `--mat`, `--pat` (splice `.mmi` or FASTA), `--out`, `--minimap2` (default `minimap2`), `--threads` (default 4).
- **Projection params:** `minimap2 -c --cs -x splice -N 50 -p 0.01 -t <threads>`; keep hits with identity ≥ 0.95 AND aligned-fraction (of the query) ≥ 0.90.
- **Assignment gate (exact):** Tier-1 SUN = assembly carries the best copy's private base at ≥1 private position → `SUN`. Tier-2 = no private position, identity margin over runner-up ≥ 0.002 → `align_fallback`. Else / Tier-3 / private-not-confirmed → `UNRESOLVED`. Single-copy family → `single_copy`.
- **Dedup:** loci with reciprocal overlap ≥ 0.50 collapse to one, keeping the highest-identity hit as `best_copy`; the second-highest overlapping identity is the `runner_up_identity`.
- **Graceful degradation:** minimap2 missing / non-zero exit on a haplotype → that haplotype yields 0 loci, WARN to stderr, continue (match `genome_projection`'s `Ok(None)` contract).
- **Output:** `<out>.parcn.tsv` (`family_id copy_id sun_tier loci_mat loci_pat parCN assign_method`) + `<out>.parcn_families.tsv` (`family_id n_copies famCN_diploid n_unresolved_loci`). `famCN_diploid = Σ parCN`.
- TDD, bite-sized, frequent commits. Focused tests only (`cargo test --lib vg_family::parcn`); the crate is large — never a full `cargo test`.

---

### Task 1: Module scaffold + `copies.fa` parser

**Files:**
- Create: `src/rustle/vg_family/parcn.rs`
- Modify: `src/rustle/vg_family/mod.rs` (register `pub mod parcn;`)

**Interfaces:**
- Produces: `pub struct Copy { pub family_id: String, pub copy_id: String, pub seq: Vec<u8> }`; `pub fn parse_copies_fa(path: &str) -> anyhow::Result<std::collections::BTreeMap<String, Vec<Copy>>>` (family_id → copies, in file order).

- [ ] **Step 1: Write the failing test** (in `parcn.rs` `mod tests`):

```rust
#[test]
fn parse_copies_fa_groups_by_family() {
    let dir = std::env::temp_dir();
    let p = dir.join(format!("parcn_copies_{}.fa", std::process::id()));
    std::fs::write(&p, ">RBMY|0|chrY:1-9|+|nexon=3\nACGTACGTA\n>RBMY|1|chrY:20-28|+|nexon=3\nACGTACGTT\n>DAZ|0|chrY:99-104|-|nexon=1\nGGGCCC\n").unwrap();
    let fams = parse_copies_fa(p.to_str().unwrap()).unwrap();
    std::fs::remove_file(&p).ok();
    assert_eq!(fams.len(), 2);
    assert_eq!(fams["RBMY"].len(), 2);
    assert_eq!(fams["RBMY"][0].copy_id, "0");
    assert_eq!(fams["RBMY"][1].seq, b"ACGTACGTT");
    assert_eq!(fams["DAZ"][0].seq, b"GGGCCC");
}
```

- [ ] **Step 2: Run, verify RED**

Run: `cargo test --lib vg_family::parcn::tests::parse_copies_fa_groups_by_family`
Expected: FAIL — module/function not found.

- [ ] **Step 3: Implement** the struct + parser + register the module:

```rust
//! Assembly-based paralog-specific copy number (parCN). OPTIONAL assembly/DNA-side supplement:
//! projects catalog copy consensuses onto phased haplotype assemblies and counts per-copy genomic loci,
//! disambiguated by deterministic SUN witnesses. Consumes only copies.fa + the assemblies; never wires
//! into the RNA-exclusive core. See docs/superpowers/specs/2026-07-14-assembly-parcn-design.md.

use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct Copy {
    pub family_id: String,
    pub copy_id: String,
    pub seq: Vec<u8>,
}

/// Parse a `gw_family_catalog` `copies.fa` (`>{family}|{copy_idx}|{chrom}:{s}-{e}|{strand}|nexon={n}`) into
/// families → copies, preserving file order. Sequence lines are concatenated and upper-cased.
pub fn parse_copies_fa(path: &str) -> anyhow::Result<BTreeMap<String, Vec<Copy>>> {
    let text = std::fs::read_to_string(path)?;
    let mut fams: BTreeMap<String, Vec<Copy>> = BTreeMap::new();
    let (mut fam, mut cid, mut seq) = (String::new(), String::new(), Vec::<u8>::new());
    let mut have = false;
    let flush = |fams: &mut BTreeMap<String, Vec<Copy>>, fam: &str, cid: &str, seq: &mut Vec<u8>| {
        if !fam.is_empty() {
            fams.entry(fam.to_string()).or_default().push(Copy { family_id: fam.to_string(), copy_id: cid.to_string(), seq: std::mem::take(seq) });
        }
    };
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if have { flush(&mut fams, &fam, &cid, &mut seq); }
            let mut it = h.split('|');
            fam = it.next().unwrap_or("").to_string();
            cid = it.next().unwrap_or("0").to_string();
            have = true;
        } else {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if have { flush(&mut fams, &fam, &cid, &mut seq); }
    Ok(fams)
}

#[cfg(test)]
mod tests {
    use super::*;
    // tests here
}
```
In `mod.rs`, add near the other `vg_family` modules: `pub mod parcn; // OPTIONAL assembly-side parCN supplement (docs/superpowers/specs/2026-07-14-assembly-parcn-design.md); never wired into the RNA-exclusive core.`

- [ ] **Step 4: Run, verify GREEN**

Run: `cargo test --lib vg_family::parcn::tests::parse_copies_fa_groups_by_family` → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/parcn.rs src/rustle/vg_family/mod.rs
git commit -m "feat(parcn): module scaffold + copies.fa parser"
```

---

### Task 2: Private positions (SUNs) from copy consensuses

**Files:**
- Modify: `src/rustle/vg_family/parcn.rs`
- Modify: `src/rustle/vg_family/copy_assign_pipeline.rs:251` (`fn banded_msa_pair` → `pub(crate) fn banded_msa_pair`)

**Interfaces:**
- Consumes: `Copy` (Task 1); `crate::vg_family::copy_assign_pipeline::banded_msa_pair(a: &[u8], b: &[u8], band: usize) -> Option<Vec<Vec<u8>>>` (returns a 2-row MSA: `[aligned_a, aligned_b]`, gap byte `b'-'`; `None` if the pair can't stay in the band).
- Produces: `#[derive(Clone,Debug,PartialEq)] pub enum Tier { T1, T2, T3, NA }`; `pub struct CopySun { pub copy_id: String, pub tier: Tier, pub private: Vec<(usize, u8)> }`; `pub fn sun_positions(copies: &[Copy], band: usize) -> Vec<CopySun>`.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn sun_positions_finds_private_snv_and_tiers() {
    // copy0 vs copy1 differ ONLY at offset 4 (A vs T) -> each has a private position there (Tier-1).
    // copy2 is identical to copy0 -> Tier-3 (indistinguishable), and offset 4 is no longer private to copy0.
    let copies = vec![
        Copy { family_id: "F".into(), copy_id: "0".into(), seq: b"ACGTAGGTCA".to_vec() },
        Copy { family_id: "F".into(), copy_id: "1".into(), seq: b"ACGTTGGTCA".to_vec() },
        Copy { family_id: "F".into(), copy_id: "2".into(), seq: b"ACGTAGGTCA".to_vec() },
    ];
    let suns = sun_positions(&copies, 8);
    let s0 = suns.iter().find(|s| s.copy_id == "0").unwrap();
    let s1 = suns.iter().find(|s| s.copy_id == "1").unwrap();
    let s2 = suns.iter().find(|s| s.copy_id == "2").unwrap();
    // copy1's 'T' at offset 4 is unique among the three -> private -> Tier-1.
    assert_eq!(s1.tier, Tier::T1);
    assert!(s1.private.iter().any(|&(p, b)| p == 4 && b == b'T'));
    // copy0 and copy2 are identical -> neither can have a private position -> Tier-3.
    assert_eq!(s0.tier, Tier::T3);
    assert_eq!(s2.tier, Tier::T3);
    assert!(s0.private.is_empty());
}

#[test]
fn sun_positions_single_copy_is_na() {
    let copies = vec![Copy { family_id: "F".into(), copy_id: "0".into(), seq: b"ACGTACGT".to_vec() }];
    let suns = sun_positions(&copies, 8);
    assert_eq!(suns[0].tier, Tier::NA);
}
```

- [ ] **Step 2: Run, verify RED**

Run: `cargo test --lib vg_family::parcn::tests::sun_positions`
Expected: FAIL — `sun_positions`/`Tier`/`CopySun` not found (and `banded_msa_pair` private).

- [ ] **Step 3: Expose the aligner + implement.** Change `copy_assign_pipeline.rs:251` `fn banded_msa_pair` to `pub(crate) fn banded_msa_pair`. Then in `parcn.rs`:

```rust
use crate::vg_family::copy_assign_pipeline::banded_msa_pair;

#[derive(Clone, Debug, PartialEq)]
pub enum Tier { T1, T2, T3, NA }

#[derive(Clone, Debug)]
pub struct CopySun { pub copy_id: String, pub tier: Tier, pub private: Vec<(usize, u8)> }

/// For copy `b` vs sibling `s`: the SET of offsets in `b` (non-gap) whose aligned base in `s` differs
/// (substitution or gap), plus `(matches, aligned_cols)` for identity. Uses the banded 2-row MSA; if the
/// pair can't be aligned in-band, treats EVERY b offset as differing (conservative: nothing private).
fn diff_offsets(b: &[u8], s: &[u8], band: usize) -> (std::collections::HashSet<usize>, usize, usize) {
    let msa = match banded_msa_pair(b, s, band) {
        Some(m) => m,
        None => return ((0..b.len()).collect(), 0, b.len().max(1)),
    };
    let (ab, asb) = (&msa[0], &msa[1]);
    let (mut boff, mut diff, mut matches, mut cols) = (0usize, std::collections::HashSet::new(), 0usize, 0usize);
    for k in 0..ab.len() {
        let (cb, cs) = (ab[k], asb[k]);
        if cb != b'-' {
            cols += 1;
            if cs == cb { matches += 1; } else { diff.insert(boff); }
            boff += 1;
        }
    }
    (diff, matches, cols.max(1))
}

/// Per-copy private positions (a SUN = an offset in copy B whose base differs from EVERY sibling) + tier.
/// Band scales with the family's copy-length spread. Threshold-free: T1 iff ≥1 private position; T3 iff a
/// sibling is ≥99.9% identical (indistinguishable); T2 otherwise; NA for a single-copy family.
pub fn sun_positions(copies: &[Copy], band: usize) -> Vec<CopySun> {
    let mut out = Vec::with_capacity(copies.len());
    for (i, b) in copies.iter().enumerate() {
        if copies.len() == 1 {
            out.push(CopySun { copy_id: b.copy_id.clone(), tier: Tier::NA, private: Vec::new() });
            continue;
        }
        // private = offsets differing from ALL siblings; max_id = closest sibling identity.
        let mut private_set: Option<std::collections::HashSet<usize>> = None;
        let mut max_id = 0.0f64;
        for (j, s) in copies.iter().enumerate() {
            if i == j { continue; }
            let (diff, matches, cols) = diff_offsets(&b.seq, &s.seq, band);
            max_id = max_id.max(matches as f64 / cols as f64);
            private_set = Some(match private_set { None => diff, Some(acc) => acc.intersection(&diff).copied().collect() });
        }
        let private_set = private_set.unwrap_or_default();
        let mut private: Vec<(usize, u8)> = private_set.iter().map(|&p| (p, b.seq[p])).collect();
        private.sort_unstable();
        let tier = if !private.is_empty() { Tier::T1 } else if max_id >= 0.999 { Tier::T3 } else { Tier::T2 };
        out.push(CopySun { copy_id: b.copy_id.clone(), tier, private });
    }
    out
}
```

- [ ] **Step 4: Run, verify GREEN**

Run: `cargo test --lib vg_family::parcn::tests::sun_positions` → PASS. Also `cargo test --lib vg_family::copy_assign_pipeline` → still PASS (visibility change only).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/parcn.rs src/rustle/vg_family/copy_assign_pipeline.rs
git commit -m "feat(parcn): per-copy private positions (SUNs) + tier from copy-consensus alignment"
```

---

### Task 3: Read assembly bases from the minimap2 `cs` tag

**Files:**
- Modify: `src/rustle/vg_family/parcn.rs`

**Interfaces:**
- Produces: `pub fn cs_bases_at(cs: &str, query: &[u8], positions: &[usize]) -> Vec<Option<u8>>` — walk a minimap2 `cs:Z:` string (query = the copy consensus that was aligned), returning the ASSEMBLY (target) base aligned to each requested query offset; `None` if that query offset is an insertion (`+`) with no target base.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn cs_bases_reads_match_and_substitution_and_insertion() {
    // query ACGTACGT (len 8). cs: 3 matches, sub (target g / query t) at q=3, 2 matches,
    // insertion of "AA" at q=6..8, tail is target-only del (does not advance query).
    // cs grammar: :N match run; *<tgt><qry> substitution; +<seq> insertion (query-only); -<seq> deletion (target-only).
    let cs = ":3*gt:2+aa-cc";
    let q = b"ACGTACGT";
    // q0 match -> assembly base 'A'(=query); q3 substitution -> assembly base 'G'(target, upper); q5 match 'C'; q6 insertion -> None.
    let got = cs_bases_at(cs, q, &[0, 3, 5, 6]);
    assert_eq!(got, vec![Some(b'A'), Some(b'G'), Some(b'C'), None]);
}
```

- [ ] **Step 2: Run, verify RED**

Run: `cargo test --lib vg_family::parcn::tests::cs_bases_reads_match`
Expected: FAIL — `cs_bases_at` not found.

- [ ] **Step 3: Implement:**

```rust
/// Walk a minimap2 `cs:Z:` short string and return, per requested QUERY offset, the aligned TARGET
/// (assembly) base. cs ops: `:N` = N matches (target base == query base); `*xy` = substitution, x = target
/// base, y = query base (advance query 1); `+seq` = insertion, query-only (target base None); `-seq` =
/// deletion, target-only (no query advance); `~` splice = target-only. Bases are returned upper-cased.
pub fn cs_bases_at(cs: &str, query: &[u8], positions: &[usize]) -> Vec<Option<u8>> {
    let want: std::collections::HashSet<usize> = positions.iter().copied().collect();
    let mut base_at: std::collections::HashMap<usize, Option<u8>> = std::collections::HashMap::new();
    let bytes = cs.as_bytes();
    let (mut k, mut qoff) = (0usize, 0usize);
    while k < bytes.len() {
        match bytes[k] {
            b':' => {
                let mut j = k + 1; let mut n = 0usize;
                while j < bytes.len() && bytes[j].is_ascii_digit() { n = n * 10 + (bytes[j] - b'0') as usize; j += 1; }
                for _ in 0..n {
                    if want.contains(&qoff) { base_at.insert(qoff, query.get(qoff).map(|b| b.to_ascii_uppercase())); }
                    qoff += 1;
                }
                k = j;
            }
            b'*' => {
                // *<target><query>
                let tgt = bytes.get(k + 1).map(|b| b.to_ascii_uppercase());
                if want.contains(&qoff) { base_at.insert(qoff, tgt); }
                qoff += 1;
                k += 3;
            }
            b'+' => {
                let mut j = k + 1;
                while j < bytes.len() && bytes[j].is_ascii_alphabetic() {
                    if want.contains(&qoff) { base_at.insert(qoff, None); }
                    qoff += 1; j += 1;
                }
                k = j;
            }
            b'-' | b'~' => { // target-only: skip the following letters/coords, no query advance
                let mut j = k + 1;
                while j < bytes.len() && bytes[j] != b':' && bytes[j] != b'*' && bytes[j] != b'+' && bytes[j] != b'-' && bytes[j] != b'~' { j += 1; }
                k = j;
            }
            _ => { k += 1; }
        }
    }
    positions.iter().map(|p| base_at.get(p).copied().flatten()).collect()
}
```

- [ ] **Step 4: Run, verify GREEN**

Run: `cargo test --lib vg_family::parcn::tests::cs_bases_reads_match` → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/parcn.rs
git commit -m "feat(parcn): read assembly bases at query offsets from the minimap2 cs tag"
```

---

### Task 4: Hybrid assignment

**Files:**
- Modify: `src/rustle/vg_family/parcn.rs`

**Interfaces:**
- Consumes: `CopySun`/`Tier` (Task 2), `cs_bases_at` (Task 3).
- Produces: `#[derive(Clone,Debug,PartialEq)] pub enum Method { Sun, AlignFallback, Unresolved, SingleCopy }`; `pub struct Locus { pub chrom: String, pub start: u64, pub end: u64, pub best_copy: String, pub identity: f64, pub runner_up_identity: f64, pub cs: String }`; `pub struct Assignment { pub copy_id: Option<String>, pub method: Method }`; `pub fn assign_locus(locus: &Locus, sun: &CopySun, best_copy_seq: &[u8]) -> Assignment`.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn assign_locus_hybrid_tiers() {
    let mk = |cs: &str, id: f64, ru: f64| Locus { chrom: "c".into(), start: 0, end: 9, best_copy: "0".into(), identity: id, runner_up_identity: ru, cs: cs.into() };
    let seq = b"ACGTAGGTCA"; // best copy's private base is 'A' at offset 4
    let sun_t1 = CopySun { copy_id: "0".into(), tier: Tier::T1, private: vec![(4, b'A')] };
    // cs shows a MATCH across offset 4 -> assembly carries 'A' -> SUN confirmed.
    assert_eq!(assign_locus(&mk(":10", 0.99, 0.90), &sun_t1, seq).method, Method::Sun);
    // cs shows a substitution at offset 4 (target g) -> assembly does NOT carry 'A' -> UNRESOLVED.
    assert_eq!(assign_locus(&mk(":4*ga:5", 0.99, 0.90), &sun_t1, seq).method, Method::Unresolved);
    // Tier-2 with a clear identity margin -> align_fallback.
    let sun_t2 = CopySun { copy_id: "0".into(), tier: Tier::T2, private: vec![] };
    assert_eq!(assign_locus(&mk(":10", 0.99, 0.90), &sun_t2, seq).method, Method::AlignFallback);
    // Tier-2 near-tie -> UNRESOLVED.
    assert_eq!(assign_locus(&mk(":10", 0.991, 0.990), &sun_t2, seq).method, Method::Unresolved);
    // Tier-3 -> UNRESOLVED. NA -> single_copy.
    let sun_t3 = CopySun { copy_id: "0".into(), tier: Tier::T3, private: vec![] };
    assert_eq!(assign_locus(&mk(":10", 0.99, 0.0), &sun_t3, seq).method, Method::Unresolved);
    let sun_na = CopySun { copy_id: "0".into(), tier: Tier::NA, private: vec![] };
    assert_eq!(assign_locus(&mk(":10", 0.99, 0.0), &sun_na, seq).method, Method::SingleCopy);
}
```

- [ ] **Step 2: Run, verify RED**

Run: `cargo test --lib vg_family::parcn::tests::assign_locus_hybrid_tiers`
Expected: FAIL — types/`assign_locus` not found.

- [ ] **Step 3: Implement:**

```rust
#[derive(Clone, Debug, PartialEq)]
pub enum Method { Sun, AlignFallback, Unresolved, SingleCopy }

#[derive(Clone, Debug)]
pub struct Locus {
    pub chrom: String, pub start: u64, pub end: u64,
    pub best_copy: String, pub identity: f64, pub runner_up_identity: f64, pub cs: String,
}

#[derive(Clone, Debug)]
pub struct Assignment { pub copy_id: Option<String>, pub method: Method }

const ALIGN_MARGIN: f64 = 0.002;

/// Hybrid assignment of a projected locus to its best copy. Tier-1: confirm the assembly carries the best
/// copy's private base (via cs) at ≥1 private position → deterministic SUN. Tier-2: assign to the best copy
/// iff its identity beats the runner-up by ≥ ALIGN_MARGIN (flagged fallback). Tier-3 / private-not-confirmed
/// / near-tie → UNRESOLVED. NA (single copy) → single_copy.
pub fn assign_locus(locus: &Locus, sun: &CopySun, best_copy_seq: &[u8]) -> Assignment {
    match sun.tier {
        Tier::NA => Assignment { copy_id: Some(locus.best_copy.clone()), method: Method::SingleCopy },
        Tier::T1 => {
            let positions: Vec<usize> = sun.private.iter().map(|&(p, _)| p).collect();
            let bases = cs_bases_at(&locus.cs, best_copy_seq, &positions);
            let confirmed = sun.private.iter().zip(bases.iter())
                .any(|(&(_, pb), got)| *got == Some(pb));
            if confirmed { Assignment { copy_id: Some(locus.best_copy.clone()), method: Method::Sun } }
            else { Assignment { copy_id: None, method: Method::Unresolved } }
        }
        Tier::T2 => {
            if locus.identity - locus.runner_up_identity >= ALIGN_MARGIN {
                Assignment { copy_id: Some(locus.best_copy.clone()), method: Method::AlignFallback }
            } else { Assignment { copy_id: None, method: Method::Unresolved } }
        }
        Tier::T3 => Assignment { copy_id: None, method: Method::Unresolved },
    }
}
```

- [ ] **Step 4: Run, verify GREEN** — `cargo test --lib vg_family::parcn::tests::assign_locus_hybrid_tiers` → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/parcn.rs
git commit -m "feat(parcn): hybrid SUN/align-fallback/unresolved locus assignment"
```

---

### Task 5: Locus dedup by reciprocal overlap

**Files:**
- Modify: `src/rustle/vg_family/parcn.rs`

**Interfaces:**
- Consumes: `Locus` (Task 4).
- Produces: `pub fn dedup_loci(mut loci: Vec<Locus>) -> Vec<Locus>` — collapse loci on the same chrom with reciprocal overlap ≥ 0.50 into one, keeping the highest-`identity` member as the survivor and setting its `runner_up_identity` to the next-highest overlapping member's identity.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn dedup_collapses_overlapping_keeps_best() {
    let mk = |s: u64, e: u64, cp: &str, id: f64| Locus { chrom: "c1".into(), start: s, end: e, best_copy: cp.into(), identity: id, runner_up_identity: 0.0, cs: ":1".into() };
    // Two heavily-overlapping hits (copy0 id .99, copy1 id .97) + one disjoint locus.
    let loci = vec![mk(1000, 2000, "0", 0.99), mk(1010, 1990, "1", 0.97), mk(50000, 51000, "3", 0.98)];
    let mut out = dedup_loci(loci);
    out.sort_by_key(|l| l.start);
    assert_eq!(out.len(), 2);
    assert_eq!(out[0].best_copy, "0");                 // highest identity wins
    assert!((out[0].runner_up_identity - 0.97).abs() < 1e-9); // runner-up recorded
    assert_eq!(out[1].best_copy, "3");
}
```

- [ ] **Step 2: Run, verify RED** — `cargo test --lib vg_family::parcn::tests::dedup_collapses_overlapping_keeps_best` → FAIL.

- [ ] **Step 3: Implement:**

```rust
fn recip_overlap(a: &Locus, b: &Locus) -> f64 {
    if a.chrom != b.chrom { return 0.0; }
    let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
    if hi <= lo { return 0.0; }
    let ov = (hi - lo) as f64;
    let la = (a.end - a.start).max(1) as f64;
    let lb = (b.end - b.start).max(1) as f64;
    (ov / la).min(ov / lb)
}

/// Collapse reciprocal-overlap ≥ 0.50 loci into one, keeping the highest-identity member (its best_copy)
/// and recording the next-highest overlapping identity as runner_up_identity (for the Tier-2 margin gate).
pub fn dedup_loci(mut loci: Vec<Locus>) -> Vec<Locus> {
    loci.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<Locus> = Vec::new();
    for l in loci {
        if let Some(k) = kept.iter_mut().find(|k| recip_overlap(k, &l) >= 0.50) {
            if k.runner_up_identity == 0.0 || l.identity > k.runner_up_identity {
                // k already has the higher identity (sorted desc); l is a runner-up candidate.
                k.runner_up_identity = k.runner_up_identity.max(l.identity.min(k.identity));
            }
        } else {
            kept.push(l);
        }
    }
    kept
}
```

- [ ] **Step 4: Run, verify GREEN** — PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/parcn.rs
git commit -m "feat(parcn): dedup overlapping loci, keep best copy + runner-up identity"
```

---

### Task 6: CIGAR/cs-retaining projection in `genome_projection.rs`

**Files:**
- Modify: `src/rustle/vg_family/genome_projection.rs`

**Interfaces:**
- Consumes: the existing private `run_minimap2_paf` (same module) — but it runs `-c` without `--cs`; add a sibling runner `run_minimap2_paf_cs` that adds `--cs` (short) so PAF carries a `cs:Z:` tag.
- Produces: `pub fn project_with_cs(consensuses: &[(String, Vec<u8>)], target: &str, min_identity: f64, min_cov: f64, minimap2: &str, threads: usize) -> anyhow::Result<Vec<ProjHit>>` where `pub struct ProjHit { pub qname: String, pub chrom: String, pub start: u64, pub end: u64, pub identity: f64, pub cov: f64, pub cs: String }` (one row per retained PAF hit; `qname` is the query FASTA header = `family|copy`).

- [ ] **Step 1: Write the failing test** (gated on minimap2, mirroring the existing projection tests):

```rust
#[test]
fn project_with_cs_returns_cs_tag() {
    if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
    // Build a tiny "genome" with two near-identical copies of a query, one carrying a SNV.
    let dir = std::env::temp_dir();
    let g = dir.join(format!("parcn_g_{}.fa", std::process::id()));
    // 300bp non-degenerate query; genome = query at 0, query+SNV at 500 (padded with Ns).
    let q: Vec<u8> = (0..300).map(|i| b"ACGT"[((i * 2654435761usize) >> 13) & 3]).collect();
    let mut snv = q.clone(); snv[150] ^= 0b100; // flip a base
    let pad = vec![b'A'; 200];
    let mut gseq = q.clone(); gseq.extend(&pad); gseq.extend(&snv); gseq.extend(&pad);
    std::fs::write(&g, format!(">c1\n{}\n", String::from_utf8_lossy(&gseq))).unwrap();
    let hits = project_with_cs(&[("F|0".into(), q.clone())], g.to_str().unwrap(), 0.90, 0.80, "minimap2", 2).unwrap();
    std::fs::remove_file(&g).ok();
    assert!(hits.len() >= 2, "should find both genomic copies");
    assert!(hits.iter().all(|h| !h.cs.is_empty()), "each hit carries a cs tag");
    assert!(hits.iter().all(|h| h.qname == "F|0"));
}
```

- [ ] **Step 2: Run, verify RED** — `cargo test --lib vg_family::genome_projection::tests::project_with_cs_returns_cs_tag` → FAIL (fn not found).

- [ ] **Step 3: Implement.** Add `run_minimap2_paf_cs` (copy `run_minimap2_paf` but insert `"--cs"` into the arg list) and `project_with_cs` + `ProjHit`. Parse each PAF line: fields 0 (qname), 5 (target=chrom), 7/8 (target start/end), 9 (residue matches), 10 (aln block len); identity = matches / block_len; cov = (qend-qstart)/qlen (fields 1,2,3); find the `cs:Z:` tag in the trailing fields. Keep hits with `identity >= min_identity && cov >= min_cov`. Return graceful empty (`Ok(vec![])`) if minimap2 fails (match `run_minimap2_paf`'s `Ok(None)` → empty). Reuse the existing `TempFile` RAII + multi-record query-FASTA writing from `project_families_batch`.

```rust
#[derive(Clone, Debug)]
pub struct ProjHit { pub qname: String, pub chrom: String, pub start: u64, pub end: u64, pub identity: f64, pub cov: f64, pub cs: String }

fn run_minimap2_paf_cs(query_path: &std::path::Path, target: &str, minimap2: &str, threads: usize) -> Result<Option<String>> {
    let out = std::process::Command::new(minimap2)
        .args(["-c", "--cs", "-x", "splice", "-N", "50", "-p", "0.01", "-t"]).arg(threads.to_string())
        .arg(target).arg(query_path).output()
        .map_err(|e| anyhow::anyhow!("minimap2 ('{minimap2}') cs projection failed: {e}"))?;
    if !out.status.success() { return Ok(None); }
    Ok(Some(String::from_utf8_lossy(&out.stdout).into_owned()))
}

pub fn project_with_cs(consensuses: &[(String, Vec<u8>)], target: &str, min_identity: f64, min_cov: f64, minimap2: &str, threads: usize) -> Result<Vec<ProjHit>> {
    let dir = std::env::temp_dir();
    let q = dir.join(format!("rustle_parcn_q_{}_{}.fa", std::process::id(), consensuses.len()));
    let _c = TempFile(q.clone());
    {
        let mut f = std::fs::File::create(&q)?;
        for (name, seq) in consensuses { if seq.is_empty() { continue; } writeln!(f, ">{name}")?; f.write_all(seq)?; writeln!(f)?; }
    }
    let qlen: std::collections::HashMap<&str, usize> = consensuses.iter().map(|(n, s)| (n.as_str(), s.len())).collect();
    let paf = match run_minimap2_paf_cs(&q, target, minimap2, threads)? { Some(p) => p, None => return Ok(Vec::new()) };
    let mut hits = Vec::new();
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 { continue; }
        let qname = f[0].to_string();
        let (qs, qe): (u64, u64) = (f[2].parse().unwrap_or(0), f[3].parse().unwrap_or(0));
        let (ts, te): (u64, u64) = (f[7].parse().unwrap_or(0), f[8].parse().unwrap_or(0));
        let (matches, blk): (f64, f64) = (f[9].parse().unwrap_or(0.0), f[10].parse().unwrap_or(1.0));
        let identity = if blk > 0.0 { matches / blk } else { 0.0 };
        let ql = *qlen.get(f[0]).unwrap_or(&1) as f64;
        let cov = if ql > 0.0 { (qe.saturating_sub(qs)) as f64 / ql } else { 0.0 };
        if identity < min_identity || cov < min_cov { continue; }
        let cs = f.iter().find_map(|t| t.strip_prefix("cs:Z:")).unwrap_or("").to_string();
        hits.push(ProjHit { qname, chrom: f[5].to_string(), start: ts, end: te, identity, cov, cs });
    }
    Ok(hits)
}
```

- [ ] **Step 4: Run, verify GREEN** — `cargo test --lib vg_family::genome_projection::tests::project_with_cs_returns_cs_tag` → PASS (or a clean no-op skip if minimap2 absent).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/genome_projection.rs
git commit -m "feat(parcn): cs-retaining genome projection variant (project_with_cs)"
```

---

### Task 7: Tabulation + TSV output

**Files:**
- Modify: `src/rustle/vg_family/parcn.rs`

**Interfaces:**
- Consumes: `Assignment`/`Method` (Task 4), `Tier` (Task 2).
- Produces: `pub struct ParcnRow { pub family_id: String, pub copy_id: String, pub tier: Tier, pub loci_mat: usize, pub loci_pat: usize, pub method: Method }`; `pub fn tabulate(family_id: &str, copies: &[Copy], suns: &[CopySun], mat: &[Assignment], pat: &[Assignment]) -> (Vec<ParcnRow>, usize)` (rows + `n_unresolved_loci`); `pub fn format_parcn_row(r: &ParcnRow) -> String`; `pub fn format_family_row(family_id: &str, rows: &[ParcnRow], n_unresolved: usize) -> String`.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn tabulate_counts_and_formats() {
    let copies = vec![
        Copy { family_id: "RBMY".into(), copy_id: "0".into(), seq: b"AAAA".to_vec() },
        Copy { family_id: "RBMY".into(), copy_id: "1".into(), seq: b"AAAT".to_vec() },
    ];
    let suns = vec![
        CopySun { copy_id: "0".into(), tier: Tier::T1, private: vec![(3, b'A')] },
        CopySun { copy_id: "1".into(), tier: Tier::T2, private: vec![] },
    ];
    let a = |cp: &str, m: Method| Assignment { copy_id: Some(cp.into()), method: m };
    let un = || Assignment { copy_id: None, method: Method::Unresolved };
    // mat: copy0 SUN once, one unresolved. pat: copy0 SUN once, copy1 fallback once.
    let mat = vec![a("0", Method::Sun), un()];
    let pat = vec![a("0", Method::Sun), a("1", Method::AlignFallback)];
    let (rows, n_unres) = tabulate("RBMY", &copies, &suns, &mat, &pat);
    let r0 = rows.iter().find(|r| r.copy_id == "0").unwrap();
    assert_eq!((r0.loci_mat, r0.loci_pat), (1, 1));       // parCN 2
    let r1 = rows.iter().find(|r| r.copy_id == "1").unwrap();
    assert_eq!((r1.loci_mat, r1.loci_pat), (0, 1));       // parCN 1
    assert_eq!(n_unres, 1);
    assert_eq!(format_parcn_row(r0), "RBMY\t0\tT1\t1\t1\t2\tSUN");
    assert_eq!(format_family_row("RBMY", &rows, n_unres), "RBMY\t2\t3\t1");
}
```

- [ ] **Step 2: Run, verify RED** — FAIL.

- [ ] **Step 3: Implement:**

```rust
#[derive(Clone, Debug)]
pub struct ParcnRow { pub family_id: String, pub copy_id: String, pub tier: Tier, pub loci_mat: usize, pub loci_pat: usize, pub method: Method }

fn tier_str(t: &Tier) -> &'static str { match t { Tier::T1 => "T1", Tier::T2 => "T2", Tier::T3 => "T3", Tier::NA => "NA" } }
fn method_str(m: &Method) -> &'static str { match m { Method::Sun => "SUN", Method::AlignFallback => "align_fallback", Method::Unresolved => "UNRESOLVED", Method::SingleCopy => "single_copy" } }

/// Count per-copy assigned loci (mat/pat), pick each copy's dominant assignment method, and total the
/// unresolved loci across both haplotypes. Copies with no assigned locus still get a row (parCN 0).
pub fn tabulate(family_id: &str, copies: &[Copy], suns: &[CopySun], mat: &[Assignment], pat: &[Assignment]) -> (Vec<ParcnRow>, usize) {
    use std::collections::HashMap;
    let tier_of: HashMap<&str, &Tier> = suns.iter().map(|s| (s.copy_id.as_str(), &s.tier)).collect();
    let mut mat_c: HashMap<String, usize> = HashMap::new();
    let mut pat_c: HashMap<String, usize> = HashMap::new();
    let mut method_of: HashMap<String, Method> = HashMap::new();
    let mut n_unres = 0usize;
    for (side, counts) in [(mat, &mut mat_c), (pat, &mut pat_c)] {
        for a in side {
            match &a.copy_id {
                Some(cp) => { *counts.entry(cp.clone()).or_insert(0) += 1; method_of.entry(cp.clone()).or_insert_with(|| a.method.clone()); }
                None => n_unres += 1,
            }
        }
    }
    let mut rows = Vec::with_capacity(copies.len());
    for c in copies {
        let tier = (*tier_of.get(c.copy_id.as_str()).unwrap_or(&&Tier::NA)).clone();
        let method = method_of.get(&c.copy_id).cloned().unwrap_or(Method::Unresolved);
        rows.push(ParcnRow { family_id: family_id.to_string(), copy_id: c.copy_id.clone(), tier,
            loci_mat: *mat_c.get(&c.copy_id).unwrap_or(&0), loci_pat: *pat_c.get(&c.copy_id).unwrap_or(&0), method });
    }
    (rows, n_unres)
}

pub fn format_parcn_row(r: &ParcnRow) -> String {
    format!("{}\t{}\t{}\t{}\t{}\t{}\t{}", r.family_id, r.copy_id, tier_str(&r.tier), r.loci_mat, r.loci_pat, r.loci_mat + r.loci_pat, method_str(&r.method))
}
pub fn format_family_row(family_id: &str, rows: &[ParcnRow], n_unresolved: usize) -> String {
    let famcn: usize = rows.iter().map(|r| r.loci_mat + r.loci_pat).sum();
    format!("{}\t{}\t{}\t{}", family_id, rows.len(), famcn, n_unresolved)
}
```

- [ ] **Step 4: Run, verify GREEN** — PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/parcn.rs
git commit -m "feat(parcn): per-copy parCN tabulation + TSV row formatting"
```

---

### Task 8: The `parcn` binary (orchestrator) + end-to-end integration test

**Files:**
- Create: `src/bin/parcn.rs`
- Modify: `Cargo.toml` (add `[[bin]] name = "parcn" path = "src/bin/parcn.rs"` if bins are not auto-discovered — check how existing bins are declared first; if `src/bin/*.rs` auto-discovers, no Cargo.toml change).

**Interfaces:**
- Consumes: everything from `crate::vg_family::parcn::*` and `genome_projection::project_with_cs`/`ProjHit`.
- Produces: the `parcn` binary. Orchestration: `parse_copies_fa` → `sun_positions` per family → for each haplotype, `project_with_cs` (all copies as `"family|copy"` queries) → group hits by family → build `Locus`es (best hit per dedup group; `runner_up_identity` from dedup) → `assign_locus` per locus → `tabulate` → write both TSVs.

- [ ] **Step 1: Write the failing integration test** in `src/bin/parcn.rs` `#[cfg(test)] mod tests` — an end-to-end run on two tiny FASTA "haplotypes" (gated on minimap2):

```rust
#[test]
fn parcn_end_to_end_two_haplotypes() {
    if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
    let dir = std::env::temp_dir();
    let tag = std::process::id();
    // Family F, two copies differing at one base (a SUN each). Each haplotype carries both copies' loci.
    let base: Vec<u8> = (0..400).map(|i| b"ACGT"[((i * 2654435761usize) >> 11) & 3]).collect();
    let mut c0 = base.clone(); let mut c1 = base.clone(); c1[200] = if c0[200] == b'A' { b'C' } else { b'A' };
    let copies_fa = dir.join(format!("parcn_e2e_copies_{tag}.fa"));
    std::fs::write(&copies_fa, format!(">F|0|c:1-1|+|nexon=1\n{}\n>F|1|c:1-1|+|nexon=1\n{}\n", String::from_utf8_lossy(&c0), String::from_utf8_lossy(&c1))).unwrap();
    let pad = vec![b'A'; 300];
    let hap = |name: &str| {
        let mut s = c0.clone(); s.extend(&pad); s.extend(&c1); s.extend(&pad);
        let p = dir.join(format!("parcn_e2e_{name}_{tag}.fa"));
        std::fs::write(&p, format!(">h_{name}\n{}\n", String::from_utf8_lossy(&s))).unwrap();
        p
    };
    let mat = hap("mat"); let pat = hap("pat");
    let out = dir.join(format!("parcn_e2e_out_{tag}"));
    // call the library orchestrator (factor main's body into `run(args) -> Result<()>` so it is testable)
    run(RunArgs { copies_fa: copies_fa.to_string_lossy().into(), mat: mat.to_string_lossy().into(), pat: pat.to_string_lossy().into(), out: out.to_string_lossy().into(), minimap2: "minimap2".into(), threads: 2 }).unwrap();
    let parcn = std::fs::read_to_string(format!("{}.parcn.tsv", out.to_string_lossy())).unwrap();
    // both copies resolved by SUN, each present on both haplotypes -> parCN 2 each.
    for cp in ["F\t0", "F\t1"] { assert!(parcn.lines().any(|l| l.starts_with(cp) && l.contains("\tSUN") && l.trim_end().ends_with("\t2\tSUN")), "copy {cp} parCN 2 via SUN"); }
    for p in [copies_fa, mat, pat] { std::fs::remove_file(p).ok(); }
    std::fs::remove_file(format!("{}.parcn.tsv", out.to_string_lossy())).ok();
    std::fs::remove_file(format!("{}.parcn_families.tsv", out.to_string_lossy())).ok();
}
```

- [ ] **Step 2: Run, verify RED** — `cargo test --bin parcn` → FAIL (no `run`/`RunArgs`).

- [ ] **Step 3: Implement `src/bin/parcn.rs`.** Define `RunArgs` + `run(RunArgs) -> anyhow::Result<()>` (the testable core), and a `clap`-parsed `main` that fills `RunArgs` and calls `run`. `run`:
  1. `let fams = parse_copies_fa(&a.copies_fa)?;`
  2. Build the query list once: `Vec<(String /*"fam|copy"*/, Vec<u8>)>` over all families' copies, and a `seq_of: HashMap<"fam|copy", Vec<u8>>`.
  3. For each haplotype target in `[&a.mat, &a.pat]`: `let hits = project_with_cs(&queries, target, 0.95, 0.90, &a.minimap2, a.threads)?;` group `hits` by family (split `qname` on `'|'`); within a family build `Locus{ chrom,start,end,best_copy=copy, identity, runner_up_identity:0, cs }` from each hit, then `dedup_loci`; for each deduped locus look up its `best_copy`'s `CopySun` + `seq_of["fam|copy"]` and call `assign_locus`. Collect `Vec<Assignment>` per family per haplotype.
  4. `let suns_by_fam: BTreeMap<fam, Vec<CopySun>> = fams.iter().map(|(f,c)| (f.clone(), sun_positions(c, band(c)))).collect();` where `band(copies)` = `(max_len - min_len) + 64` (covers copy-length spread).
  5. For each family, `tabulate(fam, copies, &suns, &mat_assign, &pat_assign)` → write `format_parcn_row` lines (+ header `family_id\tcopy_id\tsun_tier\tloci_mat\tloci_pat\tparCN\tassign_method`) and `format_family_row` (+ header `family_id\tn_copies\tfamCN_diploid\tn_unresolved_loci`).
  6. On a haplotype whose `project_with_cs` returns empty, WARN to stderr and treat as zero loci (already the natural result).

Wire copy-id key consistently as `format!("{fam}|{copy_id}")` everywhere.

- [ ] **Step 4: Run, verify GREEN** — `cargo test --bin parcn` → PASS (or clean skip without minimap2). `cargo build --bin parcn` → clean.

- [ ] **Step 5: Commit**

```bash
git add src/bin/parcn.rs Cargo.toml
git commit -m "feat(parcn): parcn binary orchestrator + end-to-end two-haplotype test"
```

---

### Task 9: Real-data validation + bench doc

**Files:**
- Create: `bench/PARCN_VALIDATION.md`

- [ ] **Step 1: Build the splice indexes once** (one-time setup, not committed):

```bash
minimap2 -x splice -d /mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi /home/juanfra/winloci_scratch/mGorGor1.mat.cur.20231122.fasta.gz
minimap2 -x splice -d /mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi /home/juanfra/winloci_scratch/mGorGor1.pat.cur.20231122.fasta.gz
```

- [ ] **Step 2: Run `parcn` on the gorilla catalog** (the `copies.fa` from a `gw_family_catalog --homology-primary` run on GGO_mm.bam; use an existing one under `winloci_scratch`/`winloci_data`). HEAVY/genome-wide → run via harness-tracked background, outputs to `/mnt/linuxdisk`:

```bash
target/release/parcn --copies-fa <ggo_catalog>.copies.fa \
  --mat /mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi \
  --pat /mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi \
  --out /mnt/linuxdisk/home/juanfraitu/winloci_data/ggo_parcn --threads 4
```

- [ ] **Step 3: Assert the checks** and record them: (a) **conservation** — for every family, `Σ parCN + n_unresolved_loci` equals the total distinct projected loci (recompute independently from the per-haplotype hit counts); (b) **RBMY ~2×** — `famCN_diploid(RBMY) ≈ 2 ×` its haploid catalog copy count (report the actual number; DAZ/GSTM as secondary CN-stable checks); (c) **method mix** — report the SUN / align_fallback / UNRESOLVED split and confirm UNRESOLVED ≈ the ~18% Tier-3 expectation. Any panic/❌ → STOP and report.

- [ ] **Step 4: Write `bench/PARCN_VALIDATION.md`** with the three results (conservation ✓, RBMY/DAZ/GSTM diploid famCN vs 2×haploid, method mix), the exact command, and a one-line note that the RNA-exclusive core is untouched. Commit.

```bash
git add bench/PARCN_VALIDATION.md
git commit -m "bench(parcn): real gorilla-catalog validation (conservation + RBMY ~2x + method mix)"
```

---

## Self-Review

- **Spec coverage:** Boundary/RNA-exclusive (Global Constraints + Task 1 module doc + Task 9 note) ✓; consensus-derived SUNs (T2) ✓; cs-tag base reading (T3) ✓; hybrid assignment (T4) ✓; dedup (T5) ✓; cs-retaining projection (T6) ✓; per-copy parCN mat/pat + roll-up output (T7) ✓; binary + inputs + graceful degradation (T8) ✓; validation conservation + RBMY 2× + method mix (T9) ✓. Tier-3 flagged-not-resolved is the UNRESOLVED path (T4) + `n_unresolved_loci` (T7) ✓.
- **Placeholder scan:** every code step carries complete code; commands are exact; no TBD/TODO.
- **Type consistency:** `Copy`, `CopySun`/`Tier`, `Locus`/`Assignment`/`Method`, `ParcnRow`, `ProjHit` are defined once (T1/T2/T4/T6/T7) and consumed with matching field/variant names; the copy key is uniformly `"{family}|{copy_id}"`; `banded_msa_pair` (pub(crate) in T2) and `project_with_cs` (T6) signatures match their call sites in T8.
