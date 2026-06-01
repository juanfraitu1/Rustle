# Haplotype-Phasing Read-to-Copy Assignment — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Assign multi-mapping long reads to their true paralog copy by the read's alleles at copy-divergent sites (haplotype phasing), make it the default `--vg` method, and remove the mis-assigning HMM-EM.

**Architecture:** A new `enumerate_diagnostic_sites` finds real divergent positions between copies by exon-pair synteny alignment (fixing `build_exon_fingerprints`, which currently finds 0 sites for genuinely-divergent copies). The already-repaired phasing scorer (`score_read_exon_fingerprint`, commit `9436e7e`) consumes those sites. Routing changes make phasing the default `--vg` path; then the HMM is deleted and shared infra moved to a `vg_graph` module.

**Tech Stack:** Rust (cargo), samtools, the multi-copy oracle (`bench/multi_copy_eval/`), gffcompare.

**Reference:** Design spec `docs/superpowers/specs/2026-05-31-haplotype-phasing-copy-assignment-design.md`.

---

## File Structure

**Phase A (phasing fix — land first):**
- Modify `src/rustle/vg.rs`: replace the inner per-ExonClass loop of `build_exon_fingerprints` with a call to new `enumerate_diagnostic_sites`; add that fn + helpers (`syntenic_exon_pairs`, `diff_exon_pair`, `exon_offset_to_genomic`); remove the uncommitted `RUSTLE_VG_HMM_SCORE_TRACE` block; keep the committed-this-session site-dedup collapse + adaptive gate; keep/commit the `n_sites==0` pileup-fallback.
- Modify `src/rustle/pipeline.rs`: ungate genome threading from `--vg-snp` (~10039); route default `--vg` to the phasing EM (~10645).
- Create `tests/regression/vg_phasing_sites.rs`: unit tests for `enumerate_diagnostic_sites` coordinate/strand invariants.
- Create `bench/multi_copy_eval/score_phasing_assignment.py`: real-family per-copy ratio check (fam 175 / fam 214 by coordinate).
- Modify `bench/multi_copy_eval/expectations.json` + `run_oracle.py`: add the phasing check.

**Phase B (HMM removal — gated on Phase A green):**
- Create `src/rustle/vg_graph/` (move `family_graph.rs`, `positional.rs`, shared parts of `rescue.rs`/`diagnostic.rs`).
- Delete `src/rustle/vg_hmm/profile.rs`, `src/rustle/vg_hmm/scorer.rs`; gut `vg_hmm/rescue.rs` HMM scoring; delete `run_pre_assembly_em_hmm`/`run_em_reweighting_hmm`/`fit_profiles_in_place`; delete `ExonClass.profile`/`per_copy_profiles`.
- Delete binaries `src/bin/{dump_family_graph,roundtrip_family_graph,loo_assembly,leave_one_out_rescue}.rs` + HMM tests.
- Update `crate::vg_hmm::` references in `vg.rs`, `pipeline.rs`, `path_extract.rs`, `transcript_filter.rs`.

---

## Conventions for every code task
- Build: `cargo build --release 2>&1 | tail -3` (must end "Finished").
- Default-isolation check after any `vg.rs`/`pipeline.rs` change: confirm the change is `--vg`-gated (no non-VG path touched). The full de-novo gate runs in Task A8.
- Region helper (used by validation): `samtools view -b ../GGO.bam <region...> > /tmp/<fam>.bam && samtools index /tmp/<fam>.bam`; per-chrom FASTA: `samtools faidx ../GGO.fasta <CHROM> > /tmp/<chrom>.fa && samtools faidx /tmp/<chrom>.fa`.
- All commits end with: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## PHASE A — Phasing site-finder + default routing + validation

### Task A0: Clean the working tree (fold in / remove this session's uncommitted diagnostics)

**Files:** Modify `src/rustle/vg.rs`

- [ ] **Step 1: Remove the `RUSTLE_VG_HMM_SCORE_TRACE` block** added during investigation (in `run_pre_assembly_em_hmm`'s E-step, the `eprintln!("[HMM-SCORE] ...")` gated block). It was diagnostic only.

- [ ] **Step 2: Keep the `n_sites == 0` pileup-fallback** already edited in `run_fingerprint_em` (the block that logs `"pileup-depth-prior fallback"` and runs the EM instead of `continue`). Confirm it is present.

- [ ] **Step 3: Build + confirm tree is intentional**

Run: `cargo build --release 2>&1 | tail -3` → Expected: "Finished".
Run: `git diff --stat src/rustle/vg.rs` → Expected: only the pileup-fallback delta remains (HMM-SCORE-TRACE gone).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg.rs
git commit -m "vg: run fingerprint-EM with pileup-prior when 0 diagnostic sites (drop investigation trace)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task A1: Genome threading on the default `--vg` path

**Files:** Modify `src/rustle/pipeline.rs:~10039-10048`

- [ ] **Step 1: Broaden the `vg_snp_genome` gate.** Change the condition from
`config.vg_mode && config.vg_snp && config.genome_fasta.is_some()` to
`config.vg_mode && config.genome_fasta.is_some()`, and update the log line to
`"[VG] Loading genome FASTA for read-to-copy phasing: {}"`. (Keep `--vg-snp` accepted as a no-op alias.)

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1 | tail -3` → Expected: "Finished".

- [ ] **Step 3: Verify mismatches now populate on the default path.** On the DAZ subset (`/tmp/daz.bam`, regenerate via the spec's reproduce block if absent) run `./target/release/rustle --vg --genome-fasta ../GGO.fasta -L /tmp/daz.bam -o /tmp/t.gtf 2>&1 | grep "Loading genome FASTA for read-to-copy phasing"`.
Expected: the line prints (genome threaded without `--vg-snp`).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "vg: thread genome into ingest for any --vg run with a FASTA (was gated on --vg-snp)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task A2: `enumerate_diagnostic_sites` — the synteny site-finder (TDD)

**Files:**
- Modify `src/rustle/vg.rs` (add `enumerate_diagnostic_sites` + helpers near `build_exon_fingerprints`)
- Test: `tests/regression/vg_phasing_sites.rs`

This is the core. The unit tests below ARE the contract — they pin the coordinate/strand behavior that is the highest-risk part. Build a small `FamilyGraph` fixture in-test (reuse the pattern in the existing `build_exon_fingerprints` unit tests at the bottom of `vg.rs`, which already construct `ExonClass`/`FamilyGraph` via `make_ec`).

- [ ] **Step 1: Write failing test — two co-linear copies, equal-length exon, one SNP**

```rust
// tests/regression/vg_phasing_sites.rs
use rustle::vg::enumerate_diagnostic_sites;            // make fn pub(crate)->pub or expose via a test shim
use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, NodeIdx};

fn ec(idx: usize, copies: &[(usize, &str, (u64,u64), char)]) -> ExonClass {
    ExonClass {
        idx: NodeIdx(idx),
        chrom: "chrT".into(),
        span: copies.iter().map(|(_,_,sp,_)| *sp).fold((u64::MAX,0),|a,b|(a.0.min(b.0),a.1.max(b.1))),
        strand: copies[0].3,
        per_copy_sequences: copies.iter().map(|(c,s,_,_)| (*c, s.as_bytes().to_vec())).collect(),
        per_copy_spans: copies.iter().map(|(c,_,sp,_)| (*c, *sp)).collect(),
        copy_specific: false, profile: None, per_copy_profiles: vec![], per_copy_cov: vec![],
    }
}

#[test]
fn snp_site_genomic_coords_forward_strand() {
    // copy0 exon @100..105 = "ACGTA", copy1 exon @500..505 = "ACCTA"  (differ at offset 2: G vs C)
    let fg = FamilyGraph {
        family_id: 0,
        nodes: vec![ec(0, &[(0,"ACGTA",(100,105),'+'),(1,"ACCTA",(500,505),'+')])],
        edges: vec![],
    };
    let fp = enumerate_diagnostic_sites(&fg, 2);
    assert_eq!(fp.n_sites, 1, "one divergent column");
    // copy0 site at genomic 102 carrying 'G'; copy1 site at genomic 502 carrying 'C'
    let s0 = &fp.per_copy_site_refs[0];
    let s1 = &fp.per_copy_site_refs[1];
    assert_eq!(s0.len(), 1); assert_eq!(s1.len(), 1);
    assert_eq!((s0[0].1, s0[0].2 as char), (102, 'G'));
    assert_eq!((s1[0].1, s1[0].2 as char), (502, 'C'));
}
```

- [ ] **Step 2: Run → fail** (`enumerate_diagnostic_sites` undefined).
Run: `cargo test --release --test vg_phasing_sites snp_site_genomic_coords_forward_strand 2>&1 | tail -5`
Expected: compile error / not found.

- [ ] **Step 3: Implement `enumerate_diagnostic_sites` (forward, gapless, equal-length first).**
Add to `vg.rs`. Reuse `ExonVariantSite`/`ExonFingerprints`. Skeleton:

```rust
/// Enumerate divergent sites between copies by syntenic exon-pair comparison.
/// Replaces the broken exon-relative inner loop of build_exon_fingerprints:
/// homologous exons may sit in SEPARATE ExonClasses (the k-mer merge fails at
/// 5-7% divergence), so we match exons across copies by genomic synteny here.
pub fn enumerate_diagnostic_sites(fg: &FamilyGraph, n_copies: usize) -> ExonFingerprints {
    // 1. Per copy: ordered exons (genomic span, strand, sequence). Collect from nodes
    //    where the copy has a per_copy_span + per_copy_sequence, sorted by span.start.
    let mut exons_by_copy: Vec<Vec<(u64,u64,char,&[u8])>> = vec![Vec::new(); n_copies];
    for node in &fg.nodes {
        for (cid, (s,e)) in &node.per_copy_spans {
            if *cid >= n_copies { continue; }
            if let Some((_, seq)) = node.per_copy_sequences.iter().find(|(c,_)| c==cid) {
                exons_by_copy[*cid].push((*s,*e,node.strand,seq.as_slice()));
            }
        }
    }
    for v in &mut exons_by_copy { v.sort_by_key(|(s,_,_,_)| *s); }

    let mut sites: Vec<ExonVariantSite> = Vec::new();
    let mut per_copy_refs: Vec<Vec<(usize,u64,u8)>> = (0..n_copies).map(|_| Vec::new()).collect();

    // 2. For each unordered copy pair, match syntenic exons and diff them.
    for i in 0..n_copies { for j in (i+1)..n_copies {
        for (ei, ej) in syntenic_exon_pairs(&exons_by_copy[i], &exons_by_copy[j]) {
            // ei/ej are (start,end,strand,seq). diff_exon_pair returns
            // (offset_in_i, base_i, offset_in_j, base_j) for each divergent column.
            for (oi, bi, oj, bj) in diff_exon_pair(ei, ej) {
                let pos_i = exon_offset_to_genomic(ei, oi);
                let pos_j = exon_offset_to_genomic(ej, oj);
                let site_idx = sites.len();
                per_copy_refs[i].push((site_idx, pos_i, bi));
                per_copy_refs[j].push((site_idx, pos_j, bj));
                sites.push(ExonVariantSite { copy_bases: vec![(i, bi), (j, bj)] });
            }
        }
    }}
    // 3. Sort + dedup per-copy site lists by genomic ref_pos (a site found in multiple
    //    pairs collapses to one entry per copy — reuse the dedup style from 9436e7e).
    for refs in &mut per_copy_refs {
        refs.sort_unstable_by_key(|&(_,p,_)| p);
        refs.dedup_by_key(|&mut (_,p,_)| p);
    }
    let n_sites = sites.len();
    if std::env::var_os("RUSTLE_VG_FP_SITE_TRACE").is_some() {
        eprintln!("[FP-TRACE] n_copies={} n_sites={}", n_copies, n_sites);
        for cid in 0..n_copies { eprintln!("[FP-TRACE]   copy {} sites={}", cid, per_copy_refs[cid].len()); }
    }
    ExonFingerprints { sites, per_copy_site_refs: per_copy_refs, n_copies, n_sites }
}

// Helpers (implement in subsequent steps with their own tests):
fn syntenic_exon_pairs<'a>(a: &'a [(u64,u64,char,&[u8])], b: &'a [(u64,u64,char,&[u8])])
    -> Vec<((u64,u64,char,&'a[u8]), (u64,u64,char,&'a[u8]))> { /* Step 6 */ unimplemented!() }
fn diff_exon_pair(a: (u64,u64,char,&[u8]), b: (u64,u64,char,&[u8]))
    -> Vec<(usize,u8,usize,u8)> { /* Step 4 */ unimplemented!() }
fn exon_offset_to_genomic(exon: (u64,u64,char,&[u8]), off: usize) -> u64 { /* Step 5 */ unimplemented!() }
```

For Step 3 minimal pass: implement `diff_exon_pair` for equal-length gapless (Step 4), `exon_offset_to_genomic` forward (Step 5), and `syntenic_exon_pairs` as positional 1:1 zip for now (Step 6 refines). Make `enumerate_diagnostic_sites` and the helpers reachable from the test (`pub`, and re-export from the crate root or test via `rustle::vg::...`).

- [ ] **Step 4: Implement `diff_exon_pair` (equal-length gapless) + test**

```rust
fn diff_exon_pair(a: (u64,u64,char,&[u8]), b: (u64,u64,char,&[u8])) -> Vec<(usize,u8,usize,u8)> {
    let (sa, sb) = (a.3, b.3);
    let mut out = Vec::new();
    if sa.len() == sb.len() {
        for o in 0..sa.len() {
            let (ba, bb) = (sa[o].to_ascii_uppercase(), sb[o].to_ascii_uppercase());
            if matches!(ba, b'A'|b'C'|b'G'|b'T') && matches!(bb, b'A'|b'C'|b'G'|b'T') && ba != bb {
                out.push((o, ba, o, bb));
            }
        }
    } else {
        // Step 7: minimizer-anchored banded fill for unequal lengths.
        out = diff_exon_pair_anchored(a, b);
    }
    out
}
```
Add a unit test `diff_equal_length_two_snps` asserting two divergent offsets on a hand-built pair.

- [ ] **Step 5: Implement `exon_offset_to_genomic` (strand-aware) + test**

```rust
// per_copy_sequences are on the TRANSCRIPT strand. For a '+' exon the sequence
// is genome-forward: genomic = start + off. For a '-' exon the sequence is the
// reverse-complement of the genome, so offset o maps to genomic end-1-o.
fn exon_offset_to_genomic(exon: (u64,u64,char,&[u8]), off: usize) -> u64 {
    let (start, end, strand, _) = exon;
    if strand == '-' { end.saturating_sub(1).saturating_sub(off as u64) } else { start + off as u64 }
}
```
Unit test `offset_to_genomic_minus_strand`: exon `(100,105,'-')`, off 2 → 102 (end-1-off = 104-2). off 0 → 104.

IMPORTANT: when strand is '-', the stored allele in `diff_exon_pair` is on the transcript strand; the scorer compares against `read.mismatches` which are GENOME-FORWARD query bases. So for a '-' exon, complement the base before storing. Encode this in `enumerate_diagnostic_sites` (complement `bi`/`bj` when the respective exon strand is '-') and add the invariant to the test below.

- [ ] **Step 6: Implement `syntenic_exon_pairs` (positional + minimizer validation) + test**

Start positional (zip by order). Then guard each pairing with a minimizer-overlap check (reuse `crate::vg_hmm::family_graph::minimizers(seq,15,10)`): if the Jaccard of two positionally-aligned exons' minimizers is below a small floor (e.g. 0.05), they are not homologs — skip (do not emit sites) and try to re-sync by searching nearby exons. For inverted copies (the two copies' exons in opposite genomic order), detect via the copies' node `strand` disagreement and reverse one list before zipping.
Unit test `inverted_copy_pairs_match`: copy0 exons ascending '+', copy1 exons the reverse-complement '-' in descending genomic order → pairs line up and the SNP maps to correct genomic coords on both.

- [ ] **Step 7: Implement `diff_exon_pair_anchored` (unequal length, banded) + test**

Minimizer-anchor the two sequences, chain co-linear anchors, and within each anchor-to-anchor block compare gaplessly the equal-length sub-run (emit SNP columns); skip indel gaps (length mismatch between anchors) — never emit a site inside an indel. Unit test `unequal_length_indel_skips_gap`: a 1bp insertion in one copy → the downstream SNP still maps to the correct genomic position (anchors absorb the shift), and no spurious sites in the gap.

- [ ] **Step 8: Run all unit tests → pass**
Run: `cargo test --release --test vg_phasing_sites 2>&1 | tail -8`
Expected: all PASS (forward SNP coords, minus-strand coords+complement, two-SNP, inverted-pair, indel-gap).

- [ ] **Step 9: Commit**

```bash
git add src/rustle/vg.rs tests/regression/vg_phasing_sites.rs
git commit -m "vg: enumerate_diagnostic_sites — find divergent sites by syntenic exon-pair alignment

Fixes diagnostic-site identification (build_exon_fingerprints found 0 sites for
genuinely-divergent copies because homologous exons land in separate ExonClasses).
Matches exons across copies by synteny+minimizers, diffs base-by-base, emits sites
in genome-forward coords with strand-correct alleles. Unit-tested for +/- strand,
inversion, and indel coordinate invariants.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task A3: Wire `enumerate_diagnostic_sites` into the EM

**Files:** Modify `src/rustle/vg.rs`

- [ ] **Step 1: Replace the site-finding call.** In `run_fingerprint_em`, change `let fp = build_exon_fingerprints(fg, n_copies);` to `let fp = enumerate_diagnostic_sites(fg, n_copies);`. Keep the `n_sites == 0` pileup-fallback and the adaptive gate untouched. (Leave `build_exon_fingerprints` in place for now; its own unit tests still pass — it becomes dead after Task A3 and is removed in Task A4 Step 1.)

- [ ] **Step 2: Build + the synthetic oracle still 100%**

Run: `cargo build --release 2>&1 | tail -3` → "Finished".
Run: `rm -f /tmp/synth_assign.attr.tsv; python3 bench/multi_copy_eval/run_oracle.py --fast --check 2>&1 | tail -3`
Expected: `ALL OBJECTIVES PASS` (synthetic Obj-4 still 100% — the synthetic copies are at separate loci with planted SNPs; the synteny finder must reproduce them).

- [ ] **Step 3: fam 175 now finds sites + flips to B>A**

Setup (if absent): `samtools view -b ../GGO.bam NC_073234.2:62950000-62965000 NC_073234.2:67144000-67162000 > /tmp/f175.bam && samtools index /tmp/f175.bam`; `samtools faidx ../GGO.fasta NC_073234.2 > /tmp/chr234.fa && samtools faidx /tmp/chr234.fa`.
Run:
```bash
RUSTLE_VG_FP_SITE_TRACE=1 ./target/release/rustle --vg --vg-no-hmm --genome-fasta /tmp/chr234.fa -L /tmp/f175.bam -o /tmp/f175.gtf 2>&1 | grep "FP-TRACE.*n_sites"
covA=$(awk -F'\t' '$3=="transcript"&&$4>=62950000&&$5<=62965000{match($9,/cov "([0-9.]+)"/,a);s+=a[1]}END{print s}' /tmp/f175.gtf)
covB=$(awk -F'\t' '$3=="transcript"&&$4>=67144000&&$5<=67162000{match($9,/cov "([0-9.]+)"/,a);s+=a[1]}END{print s}' /tmp/f175.gtf)
echo "covA=$covA covB=$covB"
```
Expected: `n_sites` now in the tens (was 0); **covB > covA** (truth direction; was inverted).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg.rs
git commit -m "vg: run_fingerprint_em uses enumerate_diagnostic_sites (real phasing signal)

fam175: 0 -> tens of diagnostic sites; per-copy ratio flips to correct B>A.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task A4: Make phasing the default `--vg` EM path

**Files:** Modify `src/rustle/pipeline.rs:~10591,10645`; Modify `src/rustle/vg.rs` (remove dead `build_exon_fingerprints` once unreferenced)

- [ ] **Step 1: Route default `--vg` to the phasing EM.** In the VG dispatch, change so that when `build_graph` (genome+graph available), `run_fingerprint_em` runs regardless of `vg_no_hmm` (i.e. drop the `do_hmm` branch as the selector; HMM call is removed in Phase B, but here we stop *selecting* it). Concretely: replace `let do_hmm = hmm_ok && !config.vg_no_hmm;` usage at the dispatch so the phasing EM is chosen for default `--vg`. Leave `run_pre_assembly_em_hmm` defined (deleted in Phase B) but no longer called.

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1 | tail -3` → "Finished" (expect dead-code warnings for HMM fns — acceptable until Phase B).

- [ ] **Step 3: Default `--vg` (no `--vg-no-hmm`) now phases fam 175 correctly**

Run the Task A3 Step-3 command but with plain `--vg` (drop `--vg-no-hmm`). Expected: covB > covA (default path now uses phasing).

- [ ] **Step 4: DAZ — no DAZ3 fabrication on the default path**

Run: `./target/release/rustle --vg --genome-fasta ../GGO.fasta -L /tmp/daz.bam -o /tmp/daz.gtf 2>/dev/null`; then
`awk -F'\t' '$3=="transcript"&&$7=="+"&&$4>=42879000&&$5<=42946000{match($9,/cov "([0-9.]+)"/,a);s+=a[1]}END{print "DAZ3 cov="s}' /tmp/daz.gtf`.
Expected: DAZ3 cov collapses toward ~0 (phasing keeps the reads at DAZ1; no fabricated copy). Record the value; encode as the gate in Task A7.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "vg: phasing-EM is the default --vg read-assignment method (HMM no longer selected)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task A5: New oracle check — real-family phasing assignment

**Files:** Create `bench/multi_copy_eval/score_phasing_assignment.py`; Modify `bench/multi_copy_eval/expectations.json`, `bench/multi_copy_eval/run_oracle.py`

- [ ] **Step 1: Write `score_phasing_assignment.py`.** Inputs: hardcoded fixture-by-coordinate panel — fam 175 (`NC_073234.2`, copyA 62950000-62965000, copyB 67144000-67162000, expect `B>A`), fam 214 (`NC_073224.2`, copyA 101546305-101591067, copyB 101593961-101625888, expect `B>A`), DAZ (`NC_073248.2`, expect DAZ3 cov < 5). For each: build the region BAM + per-chrom FASTA if absent, run `./target/release/rustle --vg --genome-fasta <chrom.fa> -L <bam> -o <gtf>`, sum per-copy `cov`, assert the expected relation. Emit JSON `{fam175:{covA,covB,pass}, fam214:{...}, daz:{daz3_cov,pass}, all_pass}`. Skip a fixture (don't fail) if `../GGO.bam`/`../GGO.fasta` absent (CI-safe).

- [ ] **Step 2: Run it against the current (post-A4) binary**

Run: `python3 bench/multi_copy_eval/score_phasing_assignment.py`
Expected: `all_pass: true` (fam175 B>A, fam214 B>A, DAZ3 cov < 5). Capture the actual numbers for expectations.

- [ ] **Step 3: Add to `expectations.json` + wire into `run_oracle.py --full --check`.** Add keys `phasing_fam175_BgtA: true`, `phasing_fam214_BgtA: true`, `phasing_daz3_cov_max: 5.0`. In `run_oracle.py`, call `score_phasing_assignment` under `--full`, compare to expectations, fail `--check` on regression (SKIPPED when GGO data absent, per the oracle's existing pattern).

- [ ] **Step 4: Determinism — run the phasing scorer twice → identical.** Run `score_phasing_assignment.py` twice; assert identical JSON.

- [ ] **Step 5: Commit**

```bash
git add bench/multi_copy_eval/score_phasing_assignment.py bench/multi_copy_eval/expectations.json bench/multi_copy_eval/run_oracle.py
git commit -m "oracle: real-family phasing-assignment check (fam175/fam214 B>A, DAZ3 no fabrication)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task A6: fam 214 + fam 136 spot-validation

- [ ] **Step 1: fam 214 ratio correct.** Region run (NC_073224.2:101540000-101630000, FASTA `/tmp/chr224.fa`); confirm copy A (101546-101591k) and copy B (101594-101626k) both assembled and the per-copy cov ratio matches phased truth (B≥A). Record.
- [ ] **Step 2: fam 136 ratio correct.** NC_073235.2:31560000-31920000 (copyA 31566817-31590386, copyB 31886302-31915543), FASTA chr235. Confirm both copies, ratio sane.
- [ ] **Step 3:** No commit (validation only); if either fails, file as a finding and STOP for controller review.

### Task A7: HARD GATE — Phase A validation

- [ ] **Step 1:** Re-run the full gate and confirm ALL pass:
  - `python3 bench/multi_copy_eval/run_oracle.py --fast --check` → `ALL OBJECTIVES PASS` (synthetic Obj-4 = 100%).
  - `python3 bench/multi_copy_eval/score_phasing_assignment.py` → `all_pass: true` (fam175/214 B>A, DAZ3 cov < 5).
  - fam 136 sane (Task A6).

### Task A8: Default de-novo isolation gate

- [ ] **Step 1: Confirm the default operating point is unchanged.** Run the de-novo benchmark (the project's standard `rustle -L GGO_19.bam` → gffcompare vs `GGO_19.gtf`, the way `run_oracle.py`'s `default_headline` / the parity bench does it).
Expected: intron-chain **96.2/91.7**, transcript **95.6/90.5** — unchanged (all Phase-A changes are `--vg`-gated).
- [ ] **Step 2:** If anything moved, STOP — a non-VG path was touched; bisect the Phase-A commits.

**>>> Phase B starts ONLY if Task A7 + A8 are green. <<<**

---

## PHASE B — Remove the HMM-EM + de-conflate shared infra (gated on Phase A)

### Task B1: Create `vg_graph` module, move shared infra

**Files:** Create `src/rustle/vg_graph/mod.rs`; move `src/rustle/vg_hmm/family_graph.rs` → `src/rustle/vg_graph/family_graph.rs`, `vg_hmm/positional.rs` → `vg_graph/positional.rs`; split shared parts of `vg_hmm/diagnostic.rs` (RescueClass framework, minimap2 helpers) into `vg_graph/diagnostic.rs`.

- [ ] **Step 1:** Move `family_graph.rs` + `positional.rs` into `vg_graph/`; create `vg_graph/mod.rs` re-exporting `FamilyGraph, ExonClass, NodeIdx, JunctionEdge, build_family_graph, recover_paralog_path, minimizers, refine_by_minimizer_jaccard, CandidateLocus, scan_genome_for_all_families, RescueClass, classify_external, cigar_has_long_indel`. Register `mod vg_graph;` in the crate root; keep `mod vg_hmm;` for now.
- [ ] **Step 2:** Update imports: replace `crate::vg_hmm::family_graph::` → `crate::vg_graph::family_graph::` (and positional/RescueClass) across `vg.rs`, `pipeline.rs`, `path_extract.rs`, `transcript_filter.rs`, `types.rs`. `enumerate_diagnostic_sites` now references `crate::vg_graph::family_graph::minimizers`.
- [ ] **Step 3: Build** → `cargo build --release 2>&1 | tail -3` → "Finished".
- [ ] **Step 4: Phase-A gate still green** (A7 commands). Commit.

```bash
git commit -am "refactor(vg): move FamilyGraph/positional/RescueClass to vg_graph (shared, not HMM)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task B2: Delete the HMM machinery

**Files:** Delete `src/rustle/vg_hmm/profile.rs`, `src/rustle/vg_hmm/scorer.rs`; Modify `src/rustle/vg.rs` (delete `run_pre_assembly_em_hmm`, `run_em_reweighting_hmm`, any `forward_against_*`/`viterbi_*`/`fit_profiles_in_place` calls); Modify `vg_graph/family_graph.rs` (delete `fit_profiles_in_place`, `ExonClass.profile`, `ExonClass.per_copy_profiles`).

- [ ] **Step 1:** Delete `profile.rs`, `scorer.rs`. Remove their `mod`/`pub use` from `vg_hmm/mod.rs`.
- [ ] **Step 2:** Delete `run_pre_assembly_em_hmm` + `run_em_reweighting_hmm` from `vg.rs` and their call sites/dispatch remnants in `pipeline.rs`. Delete `fit_profiles_in_place` and the `profile`/`per_copy_profiles` fields from `ExonClass` (and all struct-literal sites — the test helpers + rescue). Remove `RUSTLE_VG_HMM_*` env reads.
- [ ] **Step 3: Build** → fix all compile errors (dead refs). Expected end: "Finished".
- [ ] **Step 4: Phase-A gate + default-isolation gate still green** (A7 + A8). Commit.

```bash
git commit -am "refactor(vg): delete HMM-EM (profile/scorer/forward-DP/run_pre_assembly_em_hmm)

The HMM mis-assigned (inverted fam175, fabricated DAZ3); phasing replaces it. ~1.5k LOC.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task B3: Re-route Obj-2 rescue scoring off the HMM

**Files:** Modify `src/rustle/vg_hmm/rescue.rs` (or move to `vg_graph/rescue.rs`)

- [ ] **Step 1:** Replace `forward_against_family`-based candidate scoring in `run_rescue_in_memory`/`run_rescue_with_bundles` with the existing k-mer hit-fraction prefilter as the score (the rescue scaffold is unvalidated either way — only the scoring primitive changes). Remove `viterbi_against_path` usage.
- [ ] **Step 2: Build** → "Finished". **Step 3:** Phase-A gate green. Commit.

```bash
git commit -am "refactor(vg): Obj-2 rescue scoring uses k-mer identity (HMM forward-DP removed)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task B4: Delete HMM binaries + tests; collapse `vg_hmm`

**Files:** Delete `src/bin/{dump_family_graph,roundtrip_family_graph,loo_assembly,leave_one_out_rescue}.rs`; delete `tests/regression/vg_hmm_{profile,scorer,em}.rs` + HMM cases in `vg_hmm_{rescue,family_graph}.rs`; remove now-empty `vg_hmm` (fold remaining rescue into `vg_graph`).

- [ ] **Step 1:** Delete the 4 binaries (also remove their `[[bin]]` entries from `Cargo.toml`).
- [ ] **Step 2:** Delete the HMM test files; keep/port any FamilyGraph tests that exercise shared infra to `vg_graph`.
- [ ] **Step 3:** If `vg_hmm/` has nothing left but `rescue.rs`, move it to `vg_graph/rescue.rs` and delete the `vg_hmm` module.
- [ ] **Step 4: Build + full test suite** → `cargo build --release 2>&1 | tail -3`; `cargo test --release 2>&1 | tail -15` → all green.
- [ ] **Step 5:** Phase-A gate (A7) + default-isolation (A8) still green. Commit.

```bash
git commit -am "chore(vg): delete HMM diagnostic binaries + tests; collapse vg_hmm into vg_graph

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task B5: Genome-wide `--vg` sanity + docs

- [ ] **Step 1: Genome-wide `--vg` does not crash and is not net-worse than baseline.** Run `./target/release/rustle --vg --genome-fasta ../GGO.fasta -L GGO_19.bam` (chr19); gffcompare vs `GGO_19.gtf`; confirm not worse than the prior `--vg` baseline (the prior was −4 matches; phasing should be ≥ that).
- [ ] **Step 2: Update docs.** `docs/experiments/DAZ_vg_verification.md` (mark DAZ3 recovery as the false positive it is, phasing keeps reads at DAZ1); `docs/VG_OBJECTIVES_AND_ROADMAP.md` (Obj-4 now phasing-based, HMM removed); regenerate `bench/multi_copy_eval/OBJECTIVES_ASSESSMENT.md` via `--write-report`.
- [ ] **Step 3: Commit.**

```bash
git commit -am "docs(vg): phasing is the default read-to-copy method; DAZ3 documented as false positive; HMM removed

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task B6: Final review
- [ ] Dispatch a final code-reviewer over the full Phase A+B diff. Confirm: default de-novo unchanged; phasing gate green; no `vg_hmm`/HMM symbol remains; `vg_graph` holds shared infra.

---

## Self-Review (against the spec)

**Spec coverage:** §3 site-finder → A2; §5.1 algorithm → A2 Steps 3–7; §5.2 scorer reuse → A3; §5.3 genome threading → A1; §5.4 default routing → A4; §5.5 ambiguous reads → A0 Step 2 (pileup-fallback retained); §5.6 HMM removal + vg_graph → B1–B4; §6 validation gate → A5/A6/A7/A8/B5; §7 risks (coordinate tests, isolation) → A2 unit tests + A8. All covered.

**Placeholder scan:** the helper bodies in A2 are split into their own implementing steps (Steps 4–7) each with a test — not placeholders but staged TDD. No "TBD"/"handle appropriately".

**Type consistency:** `enumerate_diagnostic_sites` returns the existing `ExonFingerprints`; `ExonVariantSite { copy_bases: Vec<(usize,u8)> }` and `per_copy_site_refs: Vec<Vec<(usize,u64,u8)>>` match the structs the scorer already consumes (commit 9436e7e). Helper tuple shape `(u64,u64,char,&[u8])` is consistent across `syntenic_exon_pairs`/`diff_exon_pair`/`exon_offset_to_genomic`.

**Known judgement call:** A2 provides skeletons + tests rather than fully-compiling bodies for the banded/anchored paths (Step 7) — the tests pin the contract; the implementer iterates against `cargo test`. This is deliberate for the one genuinely-novel algorithm.
