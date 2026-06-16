# Gated read-coherence layer (`--read-coherence`) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]` checkboxes.

**Goal:** Ship read-coherent path extraction as a gated, default-off layer that adds the validated
**+580 FSM / +2,784 broad** annotated transcripts StringTie misses, while filtering the ~half noise
via degradation-aware collapse + an annotation-free realness gate (canonical + RT-switch + depth).
Byte-identical when off; `RUSTLE_PRECISE`-exempt; strictly additive over the flow floor.

**Architecture:** spec `docs/superpowers/specs/2026-06-15-read-coherence-gated-layer-design.md`.
Reuses `global_flow.rs::extract_transcripts_readchain` (read-chain engine, pipeline.rs:7284),
extends `compute_5p_degrade_folds`, adds two `genome.rs` FASTA helpers + the gate.

**Tech stack:** Rust (rustle); validation via SQANTI3 + the cached `/tmp/gw` GTFs.

**Standing constraints:** never `pkill -f rustle`; whole-genome `-L` OOMs → per-chrom serial;
`RUSTLE_PRECISE=1` byte-identical to 4705ab1; commit only when the user asks; stage explicit paths
(no `git add -A`); test builds `CARGO_PROFILE_TEST_DEBUG=0`; `RAYON_NUM_THREADS=1` for determinism.

**Key signatures (verified):**
- `GenomeIndex::fetch_sequence(&self, chrom:&str, start:u64, end:u64) -> Option<Vec<u8>>` (genome.rs:113)
- Junction coords: `donor` = left-exon-end, `acceptor` = right-exon-start, **0-based half-open**
  (intron = `[donor, acceptor)`; donor dinuc = `seq[donor..donor+2]`, acceptor dinuc = `seq[acceptor-2..acceptor]`).
- `compute_5p_degrade_folds(chains:&[(Vec<(u64,u64)>,f64,u64,u64)], strand:char, end3_tol:u64) -> Vec<Option<usize>>` (global_flow.rs:629)
- `Transcript`: `exons: Vec<(u64,u64)>`, `chrom: String`, `strand: char`, `longcov: f64`, `source: Option<String>`.

---

## Task 1: `--read-coherence` flag plumbing

**Files:** `src/rustle/main.rs` (near other long-read flags), `src/rustle/types.rs`, `src/rustle/pipeline.rs`.

- [ ] **Step 1 (failing test):** in `main.rs` tests, assert `--read-coherence` parses and sets `RUSTLE_READ_COHERENCE`. Pattern-match an existing flag test.
- [ ] **Step 2:** run → fails (flag absent).
- [ ] **Step 3:** add `#[arg(long = "read-coherence", default_value_t = false)] read_coherence: bool` in `main.rs`; when true: `std::env::set_var("RUSTLE_READ_COHERENCE","1")` AND `std::env::set_var("RUSTLE_READCHAIN","1")` (enables the engine).
- [ ] **Step 4:** run → passes. `cargo build --release`. Commit (when authorized).

## Task 2: canonical-junction helper (`genome.rs`)

**Files:** `src/rustle/genome.rs` (+ its `#[cfg(test)]`).

- [ ] **Step 1 (failing test):** with a tiny synthetic GenomeIndex whose intron reads `GT..AG` (+) and another `CT..AC` (−), assert `is_canonical_junction(chrom, donor, acceptor, '+')==true` for GT-AG, `false` for a `GG..AG` intron.
- [ ] **Step 2:** run → fails (fn missing).
- [ ] **Step 3:** implement:
```rust
/// Canonical splice site: intron [donor,acceptor) begins GT & ends AG (+ strand) or
/// begins CT & ends AC (− strand, = revcomp). strand '.' accepts either orientation.
pub fn is_canonical_junction(&self, chrom: &str, donor: u64, acceptor: u64, strand: char) -> bool {
    if acceptor < donor + 4 { return false; }
    let d = match self.fetch_sequence(chrom, donor, donor + 2) { Some(s) => s, None => return false };
    let a = match self.fetch_sequence(chrom, acceptor - 2, acceptor) { Some(s) => s, None => return false };
    let up = |b: &[u8]| -> [u8;2] { [b[0].to_ascii_uppercase(), b[1].to_ascii_uppercase()] };
    let (d, a) = (up(&d), up(&a));
    let plus  = d == *b"GT" && a == *b"AG";
    let minus = d == *b"CT" && a == *b"AC";
    match strand { '+' => plus, '-' => minus, _ => plus || minus }
}
```
- [ ] **Step 4:** run → passes. Commit (when authorized).

## Task 3: RT-switch helper (`genome.rs`)

**Files:** `src/rustle/genome.rs`.

- [ ] **Step 1 (failing test):** synthetic genome where the 8 bp before `donor` equals the 8 bp before `acceptor` → assert `is_rt_switch(...)==true`; a non-repeat case → `false`.
- [ ] **Step 2:** run → fails.
- [ ] **Step 3:** implement (direct-repeat signature; `repeat_len` default 8):
```rust
/// RT template-switch signature: a direct repeat flanking the junction (the `repeat_len` bp
/// ending at the donor equals the `repeat_len` bp ending at the acceptor). Heuristic.
pub fn is_rt_switch(&self, chrom: &str, donor: u64, acceptor: u64, repeat_len: u64) -> bool {
    if donor < repeat_len || acceptor < repeat_len { return false; }
    let up = match self.fetch_sequence(chrom, donor - repeat_len, donor) { Some(s)=>s, None=>return false };
    let dn = match self.fetch_sequence(chrom, acceptor - repeat_len, acceptor) { Some(s)=>s, None=>return false };
    !up.is_empty() && up.eq_ignore_ascii_case(&dn) && !up.iter().any(|&b| b==b'N' || b==b'n')
}
```
- [ ] **Step 4:** run → passes. Commit (when authorized).

## Task 4: degradation-aware collapse (3′ + internal)

**Files:** `src/rustle/global_flow.rs` (`compute_5p_degrade_folds` → generalize; its caller ~line 867).

- [ ] **Step 1 (failing test):** in `degrade_collapse_tests` (global_flow.rs:1877), add cases: a 3′-truncated chain (shares 5′ terminus, junctions are a prefix of a longer chain) folds into the parent; an internal-fragment chain (junctions a contiguous sub-run of a longer chain) folds; the existing 5′ case still folds; a genuinely-distinct short isoform (not a sub-path) does NOT fold.
- [ ] **Step 2:** run → fails (only 5′ folds today).
- [ ] **Step 3:** generalize: a chain B folds into super A iff B's junction list is a **contiguous sub-sequence** of A's junctions AND the shared terminus (5′ for 3′-truncation, 3′ for 5′-truncation, neither required for internal) is within `tol`. Keep behavior gated so default (5′-only) is unchanged unless `RUSTLE_READ_COHERENCE` is set (read `std::env::var_os("RUSTLE_READ_COHERENCE")`); under the flag, also fold 3′ + internal. Preserve determinism (longest-first, abundance tie-break).
- [ ] **Step 4:** run → passes; existing `degrade_collapse_tests` still green. Commit (when authorized).

## Task 5: realness gate over read-chain transcripts

**Files:** `src/rustle/pipeline.rs` (global stage, where `all_transcripts` is assembled; the genome is in scope there as the `GenomeIndex` used elsewhere).

- [ ] **Step 1 (failing test):** unit test a helper `gate_read_coherence(txs, genome, min_cov) -> Vec<Transcript>`: drops a read-chain tx (`source` starts with `"readchain"`/`"read_coherence"`) with a non-canonical junction; drops one with an RT-switch junction; drops one with `longcov < min_cov`; KEEPS a canonical, supported one; and **never touches** a `flow`/`guide` tx (returns it unchanged even if non-canonical).
- [ ] **Step 2:** run → fails (fn missing).
- [ ] **Step 3:** implement `gate_read_coherence` (in `transcript_filter.rs` or `pipeline.rs`): for each tx whose `source` marks it read-chain-derived, compute its junctions from `exons` (`(exons[i].1, exons[i+1].0)`), drop if any `!genome.is_canonical_junction(...)` OR any `genome.is_rt_switch(...)` OR `longcov < min_cov` (`RUSTLE_READ_COHERENCE_MIN_COV`, default 2). Non-read-chain tx pass through untouched. Tag survivors `source "read_coherence"`.
- [ ] **Step 4:** run → passes. Commit (when authorized).

## Task 6: wire under `--read-coherence` + invariants

**Files:** `src/rustle/pipeline.rs`.

- [ ] **Step 1:** when `RUSTLE_READ_COHERENCE` is set, after read-chain transcripts are in `all_transcripts`, call `gate_read_coherence(...)` (gated; uses the in-scope genome). Additive: read-chain tx whose intron chain duplicates a flow chain are already handled by the existing identical-transcript dedup.
- [ ] **Step 2 (invariant test):** `--read-coherence` OFF → output byte-identical to current (`rustle -G st.gtf` on chr19 unchanged); `RUSTLE_PRECISE=1` unaffected. Run, confirm.
- [ ] **Step 3:** `cargo build --release`; full lib suite green. Commit (when authorized).

## Task 7: genome-wide validation (keep the 580, drop the noise)

**Files:** `bench/` (small runner, mirror `bench/rc_gen.sh`).

- [ ] **Step 1:** per-chrom (serial, OOM-safe) `rustle -L $C.bam --read-coherence -G st_$C.gtf --genome-fasta` → `rcg_$C.gtf`; merge (re-ID by chrom).
- [ ] **Step 2:** gffcompare rcg vs RefSeq; compute **gated-rc-FSM(=) minus StringTie-FSM(=)** genome-wide (the rigorous cut, like `bench/psv_sizing/`); also broad =/c/j.
- [ ] **Step 3:** compare to the ungated `--read-chain` (10,144 raw): confirm noise dropped (raw 10,144 → ~4–5k) while the **580 FSM are retained** (≥~90%); report precision (FP) delta. Write to `bench/scorecard.md`.

## Final review
- [ ] Dispatch a final code-reviewer (spec compliance + quality + invariants: default byte-identical, additive, `RUSTLE_PRECISE` intact).
- [ ] Use superpowers:finishing-a-development-branch.
