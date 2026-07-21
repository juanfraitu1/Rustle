# Identifiability-Respecting Locus Merge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the locus merge respect read-separability — two co-located same-strand copies merge only when no read distinguishes them (χ(H)) — so the 28 "distinguishable-but-merged" Soto copies are kept separate, without touching the 36 true-K=0 merges.

**Architecture:** Compute a per-copy "has distinguishing reads" signal upstream where read `placements` already exist, store it on the `DenovoTranscript`, and have `distinct_locus_reps` consult it instead of collapsing on coordinates alone. A pure `reads_distinguish` oracle encapsulates the decision. Because neither merge call site currently has reads in scope (`reads` is dropped before `:2265`; `refine_families_exon_sum` takes no reads), the signal is carried on the copy, not threaded through the merge's callers.

**Tech Stack:** Rust (`src/rustle/vg_family/`), `cargo`, Python3 (Soto eval + floor decompose), real gorilla data in `/home/juanfra/winloci_scratch/` + `/mnt/linuxdisk/`.

## Global Constraints

- **Result-changing by design — validated on Soto, NOT byte-identity.** Success = Soto member sensitivity **> 76.2%** (up from 276/362), recovery precision **= 100%** (held), the **36 true-K=0 members still merged**, `cargo test` green.
- **No new arbitrary constant.** "Distinguishable" reuses the existing floors: conflict `min_reads` (currently 3) for unique mappers; `PSV_MIN_ALLELE_READS` (2) for a read-supported PSV; the existing `copy_junction_support` for junctions.
- **Preserve the antisense-minority rule** (`ANTISENSE_MINORITY_DENOM`) for opposite-strand pairs exactly — the change is scoped to the same-strand collapse branch.
- **Scope:** the 28 "distinguishable-but-MERGED" cases. The 10 "unseeded" (detection gap) are OUT of scope.
- **Crash rule (WSL2):** any `copy_assign`/`gw_family_catalog` run is FOREGROUND, serial, region-restricted, output to `/home/juanfra/winloci_scratch`. NO nohup/background `&`/waiter loops/pkill.
- Commit on the current branch (`dna-family-fallback`); messages end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## File Structure

- `bench/soto/merge_baseline.txt` — frozen pre-change Soto numbers (recovery/precision/K=0 families) for before/after (generated).
- `src/rustle/vg_family/read_conflict.rs` — the pure `reads_distinguish` oracle + its unit tests (co-located with the other read-conflict predicates it mirrors).
- `src/rustle/vg_family/denovo_pipeline.rs` — the `distinguishing: bool` field on `DenovoTranscript`, its population upstream (where `placements` exist), and the guard in `distinct_locus_reps`.
- `bench/mechanism/consolidation_divergences.md` — B2.3 resolution note.

---

## Task 1: Freeze the Soto baseline

**Files:**
- Create: `bench/soto/merge_baseline.txt`

- [ ] **Step 1: Capture the current Soto recovery + floor decomposition**

Run (mount must be active; `ls /mnt/linuxdisk/home/juanfraitu/winloci_data/ >/dev/null` should succeed):
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
{ echo "# merge baseline @ $(git rev-parse --short HEAD)";
  echo "## detection eval"; python3 bench/soto/soto_detection_eval.py 2>&1 | grep -E "member sensitivity|recovery precision|families fully";
  echo "## floor decomposition"; python3 bench/soto/soto_floor_decompose.py 2>&1 | grep -E "MERGED|expressed-K=0|distinguishable|silent|information floor";
} > bench/soto/merge_baseline.txt
cat bench/soto/merge_baseline.txt
```
Expected: records `member sensitivity: 276/362 = 76.2%`, `recovery precision: 100%`, and the floor buckets (28 MERGED, 36 info-floor, etc.). This is the before-picture the change is measured against.

- [ ] **Step 2: Record the specific K=0 families that must NOT change**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
echo "## K=0 families that must stay merged (must be unchanged after the fix):" >> bench/soto/merge_baseline.txt
awk -F'\t' '$NF ~ /expressed-K=0/ {print "  "$1" "$2" "$3":"$4"-"$5}' bench/soto/soto_floor_decomposition.tsv | head -20 >> bench/soto/merge_baseline.txt
tail -25 bench/soto/merge_baseline.txt
```
Expected: a list of the K=0 (info-floor) families — these have NO distinguishing reads, so the fix must leave them merged. Task 4 checks these are unchanged.

- [ ] **Step 3: Commit**

```bash
git add bench/soto/merge_baseline.txt
git commit -m "test(merge): freeze Soto baseline before the identifiability-merge change

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: The `reads_distinguish` oracle (pure function, TDD)

**Files:**
- Modify: `src/rustle/vg_family/read_conflict.rs` (add the function + `#[cfg(test)]` cases)

**Interfaces:**
- Produces: `pub fn reads_distinguish(uniq_i: usize, uniq_j: usize, shared_psv_or_junction: bool, min_reads: usize) -> bool` — returns true (keep separate) iff either copy has ≥ `min_reads` unique-mapper reads the other lacks, OR a read-supported PSV/junction separates them. The caller supplies the per-copy unique-mapper counts and the PSV/junction flag from evidence already computed upstream; this function is the pure decision so it is unit-testable in isolation.

- [ ] **Step 1: Write the failing tests**

```rust
// in src/rustle/vg_family/read_conflict.rs, inside #[cfg(test)] mod tests
#[test]
fn reads_distinguish_keeps_separate_when_unique_mappers_present() {
    // one copy has 40 unique reads (the ID_26 case): distinguishable -> keep separate
    assert!(reads_distinguish(40, 0, false, 3));
    // both sides above the floor: distinguishable
    assert!(reads_distinguish(11, 8, false, 3));
    // a read-supported PSV/junction separates them even with no unique mappers
    assert!(reads_distinguish(0, 0, true, 3));
}

#[test]
fn reads_distinguish_merges_true_k0() {
    // no unique mappers either side, no distinguishing PSV/junction -> K=0 -> merge
    assert!(!reads_distinguish(0, 0, false, 3));
    // unique support below the floor is noise, not a distinction -> merge
    assert!(!reads_distinguish(2, 1, false, 3));
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --lib reads_distinguish 2>&1 | tail -5`
Expected: FAIL (function not defined).

- [ ] **Step 3: Implement the oracle**

```rust
// src/rustle/vg_family/read_conflict.rs
/// True ⟹ the two co-located copies are DISTINGUISHABLE by reads and must be kept separate;
/// false ⟹ no read separates them (true K=0) and they may collapse. This is the χ(H) edge
/// predicate restricted to a co-located pair. No new threshold: `min_reads` is the conflict
/// floor, and the PSV/junction flag is gated upstream at `PSV_MIN_ALLELE_READS`.
pub fn reads_distinguish(uniq_i: usize, uniq_j: usize, shared_psv_or_junction: bool, min_reads: usize) -> bool {
    uniq_i >= min_reads || uniq_j >= min_reads || shared_psv_or_junction
}
```

- [ ] **Step 4: Run to verify pass**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --lib reads_distinguish 2>&1 | tail -5`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/read_conflict.rs
git commit -m "feat(merge): reads_distinguish oracle (pairwise χ(H) separability, TDD)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Carry per-copy unique-mapper evidence to the merge + apply the guard

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (the `DenovoTranscript` copy struct, its population, and `distinct_locus_reps` at `:2881` + call sites `:2265`, `:2522`)

**Interfaces:**
- Consumes: `reads_distinguish` (Task 2); the existing `family_mapq0_support` / `build_read_placements` / per-copy `uniq` machinery in `read_conflict.rs` and the family-build path.
- Produces: `distinct_locus_reps` merges a co-located same-strand pair only when `reads_distinguish` is false.

**Integration decision procedure (the wiring is the work — follow this, guarded by tests):**

The two merge call sites have NO reads in scope (`reads` is dropped before `:2265` inside `detect_homology_catalog_genome_wide`; `refine_families_exon_sum` at `:2522` takes no reads). So DO NOT thread reads into the merge. Instead carry the signal ON the copy:

1. Determine where each copy's per-copy unique-mapper count is available. Check first whether the copy transcript already carries one: `grep -n "uniq" src/rustle/vg_family/denovo_pipeline.rs` — if the copies reaching `:2265`/`:2522` already have a populated per-copy unique-mapper field, USE IT directly (no new field needed).
2. If not, add a field `pub distinguishing_uniq: usize` to the copy struct, default 0, and POPULATE it at the point where `placements` + `family_mapq0_support` are in scope (near `denovo_pipeline.rs:1426–1441`, and in `detect_homology_catalog_genome_wide` BEFORE `drop(reads)` — build placements for the reps there and record each rep's unique-mapper count). This is the only structural change; it is purely additive (a new field carried through).
3. Change `distinct_locus_reps(copies)` → `distinct_locus_reps(copies)` unchanged in signature (the evidence rides on `copies[i]`), and replace the same-strand `collapse = true` with:

```rust
// was: let collapse = if a.strand == b.strand { true } else { <antisense minority rule> };
let collapse = if a.strand == b.strand {
    // χ(H): collapse only when NO read distinguishes them (true K=0). Distinguishable copies stay separate.
    !crate::vg_family::read_conflict::reads_distinguish(
        a.distinguishing_uniq, b.distinguishing_uniq, /*shared_psv_or_junction=*/false, cfg_min_reads)
} else {
    let (lo, hi) = (a.n_reads.min(b.n_reads), a.n_reads.max(b.n_reads));
    lo.saturating_mul(ANTISENSE_MINORITY_DENOM) < hi
};
```
where `cfg_min_reads` is the conflict `min_reads` floor (thread it in as a plain `usize` param to `distinct_locus_reps`, sourced from `cfg.conflict.min_reads` — a value, not reads). Wire the PSV/junction leg (`shared_psv_or_junction`) only if that per-copy evidence is already on the copy at the merge point; if it is not, leave it `false` and note in the report that the unique-mapper leg alone is active (it covers the 28 flagged cases).

- [ ] **Step 1: Preserve the existing merge unit tests, add the new distinguishing case**

The tests at `denovo_pipeline.rs:4081` (`distinct_locus_reps_collapses_same_locus_keeps_disjoint`) and `:4095`/`:4106` (antisense) must still pass. Add a new test: two co-located same-strand copies where one has `distinguishing_uniq = 40` are KEPT SEPARATE (not collapsed), and two with `distinguishing_uniq = 0` still collapse. Write it against whatever field/param the decision procedure lands on.

```rust
#[test]
fn distinct_locus_reps_keeps_distinguishable_colocated_copies_separate() {
    // two overlapping same-strand copies; one carries 40 unique-mapper reads -> must stay 2 loci
    let a = /* DenovoTranscript at chrX:100-200, strand +, distinguishing_uniq = 40 */;
    let b = /* DenovoTranscript at chrX:150-250, strand +, distinguishing_uniq = 0 */;
    let loci = distinct_locus_reps(vec![a, b], /*min_reads=*/3);
    assert_eq!(loci.len(), 2, "a copy with 40 unique reads must not be merged");
}
#[test]
fn distinct_locus_reps_still_merges_indistinguishable_colocated_copies() {
    // two overlapping same-strand copies, neither with unique support -> still 1 locus (K=0)
    let a = /* chrX:100-200, +, distinguishing_uniq = 0 */;
    let b = /* chrX:150-250, +, distinguishing_uniq = 0 */;
    let loci = distinct_locus_reps(vec![a, b], 3);
    assert_eq!(loci.len(), 1, "indistinguishable co-located copies still merge (K=0)");
}
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --lib distinct_locus_reps 2>&1 | tail -8`
Expected: the two new tests FAIL (guard not yet applied / signature not updated); the 3 existing ones still describe current behavior.

- [ ] **Step 3: Apply the decision procedure above** — add/populate the per-copy evidence field, update `distinct_locus_reps`'s signature to take `min_reads: usize`, replace the same-strand collapse with the `reads_distinguish` guard, and update the two call sites (`:2265`, `:2522`) + the 3 `distinct_locus_reps`/`refine_families_exon_sum` call chains to pass `cfg.conflict.min_reads`.

- [ ] **Step 4: Build + run all merge tests**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -3
cargo test --lib distinct_locus_reps 2>&1 | tail -8
cargo test --lib 2>&1 | tail -4
```
Expected: builds; the 2 new tests PASS; the 3 existing merge tests still PASS (disjoint kept separate, antisense preserved); full suite green except the 1 known pre-existing deleted-fixture failure.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(merge): distinct_locus_reps merges co-located copies only when reads can't distinguish them (χ(H))

Same-strand co-located copies now collapse iff no read distinguishes them (unique-mapper
floor), replacing the unconditional coordinate collapse. Per-copy unique-mapper evidence is
carried on the copy from where placements exist; antisense rule preserved.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Validate on Soto (recovery up, precision held, K=0 untouched)

**Files:** none created; a verification-only task producing an after-picture.

- [ ] **Step 1: Regenerate the Soto detection catalog with the changed binary**

The Soto eval reads precomputed detected-family TSVs; regenerate them with the new binary first. Run the catalog build the eval consumes (foreground, crash rule) — use the same command that produced `bench/soto/a119b_detected_families.tsv` originally (check `bench/soto/README` or the eval's header for the exact invocation; it is a `gw_family_catalog`/`copy_assign` run on `A119b.t2t.bam` vs `chm13v2.0.fa`). If the eval consumes committed TSVs that must be regenerated, regenerate them into `bench/soto/` with the new binary.

- [ ] **Step 2: Re-run the eval + floor decomposition, compare to baseline**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
echo "=== AFTER ==="; python3 bench/soto/soto_detection_eval.py 2>&1 | grep -E "member sensitivity|recovery precision|families fully"
python3 bench/soto/soto_floor_decompose.py 2>&1 | grep -E "MERGED|distinguishable|information floor"
echo "=== BASELINE (before) ==="; grep -E "member sensitivity|recovery precision|MERGED" bench/soto/merge_baseline.txt
```
Expected — the acceptance gate: member sensitivity **> 76.2%** (up), recovery precision **= 100%** (held), the "distinguishable-but-MERGED" bucket **shrunk**. If precision DROPPED below 100%, the guard over-split — STOP and report (the min_reads floor or the PSV leg needs adjustment); do not accept a precision regression.

- [ ] **Step 3: Confirm the 36 K=0 families are unchanged**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# every K=0 family recorded in merge_baseline.txt must STILL be a single merged locus (not split)
grep "must stay merged" -A20 bench/soto/merge_baseline.txt | grep "ID_" | while read _ id gene loc; do
  cnt=$(awk -F'\t' -v id="$id" '$1==id' bench/soto/a119b_detected_copies.tsv | wc -l)
  echo "$id $gene: $cnt detected copies (baseline expectation: unchanged/merged)"
done | head -20
```
Expected: the K=0 families' copy counts are unchanged vs baseline (they have no distinguishing reads, so the guard leaves them merged). Any K=0 family that gained a copy is a false split — investigate.

- [ ] **Step 4: Commit the after-picture**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 bench/soto/soto_detection_eval.py 2>&1 | grep -E "member sensitivity|recovery precision|families fully" > bench/soto/merge_after.txt
git add bench/soto/merge_after.txt bench/soto/a119b_detected_*.tsv bench/soto/soto_floor_decomposition.tsv
git commit -m "test(merge): Soto after-picture — recovery up, precision held, K=0 unchanged

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Resolve B2.3 in the divergence log + re-baseline the B gate

**Files:**
- Modify: `bench/mechanism/consolidation_divergences.md`
- Modify: `bench/mechanism/byte_identity_baseline.txt` (re-freeze to the improved output)

- [ ] **Step 1: Record that B2.3 is resolved by read-separability**

Append to `bench/mechanism/consolidation_divergences.md`:
```
## B2.3 RESOLVED (2026-07-21, identifiability-merge)
The ≥2-loci over-merge (any-overlap vs reciprocal-50%) is not an overlap-threshold choice: it is
resolved by read-separability — distinct_locus_reps now merges co-located copies only when no read
distinguishes them (χ(H); see docs/superpowers/specs/2026-07-21-identifiability-merge-design.md).
Coordinate-threshold-free. Soto member recovery: <before> → <after>.
```
(fill in the before/after numbers from `merge_baseline.txt` / `merge_after.txt`).

- [ ] **Step 2: Re-baseline the B byte-identity gate to the improved output**

The change deliberately alters the catalog, so the B gate now fails by design. Re-freeze it to the new, better output (only AFTER Task 4 confirmed recovery↑/precision-held):
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -2
bash bench/mechanism/byte_identity_gate.sh freeze
head -1 bench/mechanism/byte_identity_baseline.txt   # new @ <sha>
```
Expected: a new baseline at the current commit. Note in the commit that this re-baseline reflects the intended merge improvement, not a regression.

- [ ] **Step 3: Commit**

```bash
git add bench/mechanism/consolidation_divergences.md bench/mechanism/byte_identity_baseline.txt
git commit -m "docs(merge): resolve B2.3 via read-separability; re-baseline B gate to improved catalog

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review Notes (for the planner)

- **Spec coverage:** §1 scope (28 merged, K=0 untouched) → Tasks 1/3/4. §2 criterion/oracle → Task 2 + Task 3 guard. §3 evidence-on-the-copy → Task 3 decision procedure. §4 validation → Task 4. §6 + B2.3 → Task 5. All covered.
- **The honest hard part is Task 3** — the wiring. The plan does NOT pretend the exact threading is known: it gives a decision procedure (use existing per-copy `uniq` if present, else add a field populated where `placements` exist), the known constraints (reads dropped at `:2265`, absent in `refine_families_exon_sum`), and TWO concrete unit tests as the guardrail so the wiring is verifiable however it lands. This is the legitimate way to plan an integration whose exact shape depends on runtime tracing.
- **Placeholder scan:** the Task 3 test bodies use `/* … */` for the `DenovoTranscript` literal construction because that struct's full field set is large and constructed via existing test helpers — the implementer fills the copy fields the same way the neighbouring tests at `:4081`–`:4114` already do. The behavioral assertions (`loci.len() == 2` / `== 1`) are concrete. This is the one intentional fill-in, bounded to matching an existing test-construction pattern in the same file.
- **No new number:** the oracle takes `min_reads` from the existing conflict floor; no constant is introduced.
- **Type consistency:** `reads_distinguish(uniq_i, uniq_j, shared_psv_or_junction, min_reads) -> bool` is used identically in Task 2 (definition) and Task 3 (call). `distinct_locus_reps` gains a `min_reads: usize` param in Task 3, reflected in its new tests.
