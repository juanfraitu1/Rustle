# Deliverable B — Byte-Identical Code Consolidation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the duplicate implementations deliverable A disclosed (two refiners, two ≥2-loci predicates, three EM paths) and fix config hygiene (CLI flags, hardcoded paths) — **without changing any scientific number.**

**Architecture:** An empirical byte-identity gate (real-data md5 baseline) is built first and is the test for every change. B1 (config hygiene) is byte-identical by construction. Each B2 consolidation is *attempted*, then the gate decides: md5-identical → land it; md5 differs → revert, document the divergence as a finding, defer to deliverable C. B never merges a change that alters an output.

**Tech Stack:** Rust (existing binaries `target/release/{copy_assign,gw_family_catalog,family_define,asj}`), `cargo`, bash + `md5sum`, `samtools`/`minimap2` on PATH, real gorilla data in `/home/juanfra/winloci_scratch/`.

## Global Constraints

- **Hard invariant: B never changes a number.** Every landed commit reproduces the frozen md5 baseline. A change that alters any output md5 is NOT merged (§escalation).
- **Byte-identity gate is the test.** The gate corpus: `copy_assign` GSTM (`NC_073224.2:129160000-129230000`), `gw_family_catalog` on a fixed small region list, `asj` ABCC4 — on `/home/juanfra/winloci_scratch/GGO_mm.bam` + `GGO.fasta`.
- **Crash rule (WSL2):** every `copy_assign`/`gw_family_catalog`/`asj` run is FOREGROUND, serial, region-restricted, output to `/home/juanfra/winloci_scratch`. NO `nohup`, NO background `&`, NO waiter loops, NO `pkill`.
- **`cargo test` stays green throughout.** Consolidations are refactors; existing tests must not break.
- **Escalation on divergence:** revert the routing (keep both impls), append the specific divergence to `bench/mechanism/consolidation_divergences.md`, commit that documentation. The task is complete either way (landed OR documented).
- **Out of scope:** the 5 ε values (deliverable C), the assembler carve, CLI flags for the 137 inert constants.
- Commit messages end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Work on the current branch (`dna-family-fallback`); do not create a new branch.

---

## File Structure

- `bench/mechanism/byte_identity_gate.sh` — freeze/check the real-data md5 baseline (the gate).
- `bench/mechanism/byte_identity_baseline.txt` — frozen md5s + base commit (generated).
- `bench/mechanism/consolidation_divergences.md` — log of caught divergences (may stay near-empty).
- `src/rustle/vg_family/asj_strand_bias.rs`, `asj_verify.rs`, `driver.rs`, `o2_materialize.rs` — hardcoded-path fixes (B1a).
- `src/rustle/vg_family/copy_assign_pipeline.rs`, `o2_materialize.rs`, `src/bin/copy_assign.rs` — POA_CAP/READ_CAP flags (B1b).
- `src/bin/family_define.rs` — `--legacy` gating (B2.1).
- `src/rustle/vg_family/copy_assign_pipeline.rs` / `em_copy_assign.rs` — EM routing (B2.2).
- `src/rustle/vg_family/denovo_pipeline.rs` / `driver.rs` — ≥2-loci predicate (B2.3).
- `src/rustle/vg_family/family_definition.rs` / `family_split.rs` — refiner (B2.4).
- `docs/ONE_METHOD.md` — Consistency table update (final task).

---

## Task 1: The byte-identity gate harness

**Files:**
- Create: `bench/mechanism/byte_identity_gate.sh`
- Create (generated): `bench/mechanism/byte_identity_baseline.txt`

**Interfaces:**
- Produces: `byte_identity_gate.sh freeze` writes the baseline; `byte_identity_gate.sh check` exits 0 iff every current md5 equals the baseline, else non-zero naming the drifted file.

- [ ] **Step 1: Write the gate script**

```bash
# bench/mechanism/byte_identity_gate.sh
#!/usr/bin/env bash
# Byte-identity gate for deliverable B. FOREGROUND, serial, winloci_scratch (crash rule).
# Usage: byte_identity_gate.sh freeze   # write baseline
#        byte_identity_gate.sh check    # diff current md5s vs baseline, non-zero on mismatch
set -euo pipefail
MODE="${1:?freeze|check}"
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
BIN=$RUSTLE/target/release
SCRATCH=/home/juanfra/winloci_scratch
BAM=$SCRATCH/GGO_mm.bam
FASTA=$SCRATCH/GGO.fasta
BASELINE=$RUSTLE/bench/mechanism/byte_identity_baseline.txt
WORK=$SCRATCH/bgate
mkdir -p "$WORK"; cd "$WORK"

# --- the fixed corpus (region-restricted, foreground) ---
"$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" \
  --region NC_073224.2:129160000-129230000 --homology-primary --min-copies 2 \
  --dump-psv --out gstm >/dev/null 2>&1
"$BIN/gw_family_catalog" --bam "$BAM" --fasta "$FASTA" \
  --regions "$RUSTLE/bench/mechanism/gate_regions.txt" --out cat >/dev/null 2>&1 || true

# --- collect md5s of the scientific outputs ---
collect(){
  for f in gstm.assignments.tsv gstm.families.tsv gstm.quant.tsv gstm.psv_cols.tsv \
           cat.families.tsv cat.copies.tsv; do
    [ -f "$f" ] && printf '%s  %s\n' "$(md5sum < "$f" | cut -d' ' -f1)" "$f"
  done
}

if [ "$MODE" = freeze ]; then
  { echo "# baseline @ $(git -C "$RUSTLE" rev-parse --short HEAD)"; collect; } > "$BASELINE"
  echo "froze $(grep -vc '^#' "$BASELINE") md5s -> $BASELINE"
else
  cur=$(collect)
  base=$(grep -v '^#' "$BASELINE")
  if [ "$cur" = "$base" ]; then echo "GATE PASS: all md5 identical to baseline"; exit 0
  else echo "GATE FAIL: md5 drift vs baseline:"; diff <(echo "$base") <(echo "$cur") || true; exit 1; fi
fi
```

- [ ] **Step 2: Create the catalog region list**

```bash
# bench/mechanism/gate_regions.txt  (small fixed set that exercises the refiner + >=2-loci path)
NC_073224.2:129160000-129230000
NC_073248.2:19558214-19739408
```

- [ ] **Step 3: Freeze the baseline**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
chmod +x bench/mechanism/byte_identity_gate.sh
bash bench/mechanism/byte_identity_gate.sh freeze
cat bench/mechanism/byte_identity_baseline.txt
```
Expected: a `# baseline @ <sha>` line plus several `<md5>  <file>` lines (at least the gstm.* set). If `gw_family_catalog` produced no output, the `|| true` keeps the freeze going with the copy_assign md5s only — note it and proceed.

- [ ] **Step 4: Self-test the gate (freeze → check passes; corrupt → check fails)**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/mechanism/byte_identity_gate.sh check   # must PASS immediately after freeze
# prove it can fail: perturb one output and re-check
echo "corrupt" >> /home/juanfra/winloci_scratch/bgate/gstm.assignments.tsv
bash bench/mechanism/byte_identity_gate.sh check; echo "exit=$?"   # must FAIL, exit non-zero
# restore by re-running the corpus (check re-runs copy_assign, overwriting the corrupted file)
bash bench/mechanism/byte_identity_gate.sh check   # PASS again
```
Expected: PASS, then FAIL (exit 1) after corruption, then PASS again (the re-run regenerates the file). This proves the gate actually detects drift.

- [ ] **Step 5: Commit**

```bash
git add bench/mechanism/byte_identity_gate.sh bench/mechanism/gate_regions.txt bench/mechanism/byte_identity_baseline.txt
git commit -m "test(consolidation): byte-identity gate + frozen real-data md5 baseline

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: Seed the divergence log**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cat > bench/mechanism/consolidation_divergences.md <<'EOF'
# Consolidation divergences (deliverable B)

Consolidations attempted whose byte-identity gate FAILED — kept both impls, deferred to
deliverable C. Empty means every attempted consolidation was byte-identical.

| B2 item | canonical target | what diverged (families/rows) | deferred to |
|---|---|---|---|
EOF
git add bench/mechanism/consolidation_divergences.md
git commit -m "docs(consolidation): seed divergence log

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: B1a — remove the 6 hardcoded winloci_scratch path defaults

**Files:**
- Modify: `src/rustle/vg_family/asj_strand_bias.rs:58`, `asj_verify.rs:75`, `driver.rs:103`, `o2_materialize.rs:59,255,~421`

**Interfaces:**
- Produces: no shipped-binary code path silently falls back to a `/home/juanfra/winloci_scratch` string; callers pass paths explicitly or get a clear error.

- [ ] **Step 1: Locate every hardcoded default and its consumers**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
grep -rn '"/home/juanfra/winloci_scratch' src/rustle/vg_family/*.rs src/bin/*.rs
# for each, find who reads the const:
for c in DEFAULT_BAM FAM_TSV_DEFAULT BAM_DEFAULT; do echo "== $c =="; grep -rn "$c" src/; done
```
Expected: the 6 sites from the spec. Determine for each whether it is (a) a genuine default consumed when a CLI arg is absent, or (b) a doc-comment/test-only constant. Only (a) needs fixing; (b) can stay if it never reaches a shipped run (note which).

- [ ] **Step 2: Replace each genuine default with error-if-unset**

For a constant consumed as a fallback, change the consumer from
`let bam = args.bam.unwrap_or(DEFAULT_BAM)` (or equivalent) to require the value:
```rust
// pattern — make the path required, no silent scratch fallback
let bam = args.bam.as_deref()
    .expect("--bam is required (no hardcoded winloci_scratch fallback)");
```
If the constant is only referenced in a doc-comment or `#[cfg(test)]`, leave it but add `// test-only; not a shipped default` so a grep reviewer sees it is inert. Delete constants with zero non-test consumers.

- [ ] **Step 3: Build + test**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -5
cargo test --lib 2>&1 | tail -5
```
Expected: builds; tests green. If a test relied on the hardcoded default, update it to pass an explicit path.

- [ ] **Step 4: Gate check (must stay byte-identical)**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/mechanism/byte_identity_gate.sh check
```
Expected: `GATE PASS`. (The thesis corpus always passes explicit paths, so removing the fallback cannot change its output.) If it FAILS, a shipped path was actually depending on the hardcoded default — investigate before proceeding.

- [ ] **Step 5: Verify no silent fallback remains + commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
grep -rn '"/home/juanfra/winloci_scratch' src/ | grep -viE 'test|///|//!' && echo "REMAINING (check each)" || echo "no silent scratch defaults remain"
git add src/rustle/vg_family/asj_strand_bias.rs src/rustle/vg_family/asj_verify.rs src/rustle/vg_family/driver.rs src/rustle/vg_family/o2_materialize.rs
git commit -m "refactor(consolidation): B1a remove hardcoded winloci_scratch path defaults

Paths are now required args / error-if-unset, no silent single-machine fallback.
Byte-identity gate passes (thesis runs pass explicit paths).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: B1b — expose POA_CAP and READ_CAP as CLI flags (unchanged defaults)

**Files:**
- Modify: `src/rustle/vg_family/copy_assign_pipeline.rs:524` (POA_CAP), `src/rustle/vg_family/o2_materialize.rs:250` (READ_CAP), `src/bin/copy_assign.rs` (Args struct + threading)

**Interfaces:**
- Produces: `copy_assign --poa-cap <N>` (default 20000) and `--read-cap <N>` (default 6000); the values thread to the two use sites. Defaults equal the current constants → byte-identical.

- [ ] **Step 1: Add the flags to the Args struct**

In `src/bin/copy_assign.rs` Args struct (near the other `#[arg(long, default_value_t = ...)]` fields, e.g. after `--min-copies`):
```rust
    /// poasta→minimap2 fallback length cap (bp). Above this a copy pair uses minimap2, not poasta.
    #[arg(long, default_value_t = 20_000)]
    poa_cap: usize,

    /// Per-family multimapper read-pool cap (o2_materialize). Reads beyond this are dropped by BAM order.
    #[arg(long, default_value_t = 6_000)]
    read_cap: usize,
```

- [ ] **Step 2: Thread the values to the use sites**

`POA_CAP` at `copy_assign_pipeline.rs:524` is a `const` inside a fn — change that fn to take a `poa_cap: usize` parameter (default supplied by the caller chain from `args.poa_cap`); replace the `const POA_CAP` use with the parameter. Same for `READ_CAP` in `o2_materialize.rs` (thread `read_cap` through the fetch path). Keep the `const` as the default the CLI sets, so no call site outside the CLI changes behavior.

Minimal, behavior-preserving shape:
```rust
// copy_assign_pipeline.rs — was: const POA_CAP: usize = 20_000; ... if len > POA_CAP
pub fn ...(..., poa_cap: usize) { ... if len > poa_cap { /* minimap2 */ } ... }
```

- [ ] **Step 3: Build + test**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -5
cargo test --lib 2>&1 | tail -5
```
Expected: builds; tests green.

- [ ] **Step 4: Gate check + explicit-default equivalence**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/mechanism/byte_identity_gate.sh check
# and confirm passing the default explicitly is identical to omitting it:
cd /home/juanfra/winloci_scratch
/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign --bam GGO_mm.bam --fasta GGO.fasta \
  --region NC_073224.2:129160000-129230000 --homology-primary --min-copies 2 --poa-cap 20000 --read-cap 6000 --out gstm_explicit
md5sum gstm_explicit.assignments.tsv bgate/gstm.assignments.tsv
```
Expected: `GATE PASS`; the two md5s match (explicit default ≡ omitted).

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add src/bin/copy_assign.rs src/rustle/vg_family/copy_assign_pipeline.rs src/rustle/vg_family/o2_materialize.rs
git commit -m "feat(consolidation): B1b --poa-cap/--read-cap flags, defaults unchanged (byte-identical)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: B2.1 — family_define → --legacy parity fixture

**Files:**
- Modify: `src/bin/family_define.rs`

**Interfaces:**
- Produces: `family_define` clearly labeled as a legacy parity fixture (help text + a `--legacy` acknowledgement or a stderr notice), `gw_family_catalog` documented as the sole catalog. No behavior change to either binary → byte-identical.

- [ ] **Step 1: Confirm family_define is not on any live thesis path**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
grep -rn "family_define" bench/OBJECTIVES_FLOW.md docs/ONE_METHOD.md bench/mechanism/ 2>/dev/null | head
head -60 src/bin/family_define.rs | grep -iE "legacy|frozen|parity|oracle|python"
```
Expected: family_define reproduces a frozen Python catalog from precomputed TSVs (per ONE_METHOD). Confirm the live catalog is `gw_family_catalog`.

- [ ] **Step 2: Relabel family_define as a legacy fixture (docs/help only — no logic change)**

Add to `src/bin/family_define.rs`'s command doc-comment / `about`:
```rust
// at the #[command(...)] about text
about = "LEGACY parity fixture: reproduces the frozen Python catalog from precomputed TSVs. \
The shipped genome-wide catalog is `gw_family_catalog` (exon-sum asm20 id>=0.80 edge). \
Kept for regression parity only; not the live family definition. See rustle_mechanism.html."
```
Do NOT change any computation — this task is documentation/labeling, so it is byte-identical trivially. (Actually retiring the binary is out of scope; this just removes the "which is the real one?" ambiguity.)

- [ ] **Step 3: Build + gate check**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -3
bash bench/mechanism/byte_identity_gate.sh check
```
Expected: builds; `GATE PASS` (no computation changed).

- [ ] **Step 4: Commit**

```bash
git add src/bin/family_define.rs
git commit -m "docs(consolidation): B2.1 label family_define as legacy parity fixture (gw_family_catalog is sole catalog)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: B2.2 — attempt EM consolidation (soft_quantify_em → em_assign_family)

**Files:**
- Modify: `src/rustle/vg_family/copy_assign_pipeline.rs` (4 `soft_quantify_em` call sites: :1331, :1553, :2375, :2430)
- Reference: `src/rustle/vg_family/em_copy_assign.rs:226` `em_assign_family`

**Interfaces:**
- Consumes: `soft_quantify_em(read_obs, alleles, QUANT_ERROR=0.01, 100)` vs `em_assign_family(read_obs, copy_alleles, &[], &[], &params, 1e-6, 200)`.
- Produces: EITHER a landed consolidation (gate passed) OR a documented divergence (gate failed → reverted).

**Expectation:** the two use different error rates (0.01 vs 0.003) and different signatures/convergence — the gate will very likely FAIL on `gstm.quant.tsv`. That is the expected, correct outcome; document and defer to C.

- [ ] **Step 1: Attempt the routing on ONE call site first (the default quant path)**

Route `copy_assign_pipeline.rs:1331`'s `soft_quantify_em(&read_obs_for_em, &surviving_alleles, QUANT_ERROR, 100)` to `em_assign_family` with equivalent inputs (map the args; `em_assign_family` returns richer output — extract the abundance vector). Leave the other 3 call sites for now. Keep `soft_quantify_em` defined (do not delete).

- [ ] **Step 2: Build + test**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -5
cargo test --lib 2>&1 | tail -5
```
Expected: builds; tests green (or a quant-related test changes — note it).

- [ ] **Step 3: Gate check — this decides the outcome**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/mechanism/byte_identity_gate.sh check
```
- **If GATE PASS** (unexpected but great): the two EM paths were equivalent on this corpus. Route the other 3 call sites too, re-run the gate, and if still PASS, commit as a landed consolidation:
  ```bash
  git add src/rustle/vg_family/copy_assign_pipeline.rs
  git commit -m "refactor(consolidation): B2.2 EM paths -> em_assign_family (gate PASS, byte-identical)

  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```
- **If GATE FAIL** (expected): go to Step 4.

- [ ] **Step 4 (if gate failed): revert, document, defer**

Capture WHAT diverged, then revert:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# show the diff the gate found (its output already printed it); capture the quant rows:
diff <(sort /home/juanfra/winloci_scratch/bgate/gstm.quant.tsv) <(git show HEAD:/dev/null 2>/dev/null; echo) 2>/dev/null | head
git checkout -- src/rustle/vg_family/copy_assign_pipeline.rs   # revert the routing
cargo build --release 2>&1 | tail -3
bash bench/mechanism/byte_identity_gate.sh check                # must PASS again after revert
```
Append to `bench/mechanism/consolidation_divergences.md` a row:
```
| B2.2 EM | em_assign_family | gstm.quant.tsv abundances differ (soft_quantify_em uses QUANT_ERROR=0.01 + 100-iter; em_assign_family uses AssignParams error_rate=0.003 + eps=1e-6/200-iter — different error model AND convergence). Entangled with the 5-epsilon inconsistency. | deliverable C |
```
Commit the documentation (NOT the reverted code):
```bash
git add bench/mechanism/consolidation_divergences.md
git commit -m "docs(consolidation): B2.2 EM paths DIVERGE (different error model) -> defer to deliverable C

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
Expected: gate PASSES after revert (baseline restored); divergence recorded. Task complete (documented outcome).

---

## Task 6: B2.3 — attempt ≥2-loci predicate consolidation (distinct_locus_reps → distinct_loci)

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs:2265,2522` (call sites of `distinct_locus_reps`)
- Reference: `distinct_locus_reps` (denovo_pipeline.rs:2881, any-overlap), `distinct_loci` (driver.rs:212, reciprocal-50%)

**Interfaces:**
- Produces: EITHER landed (gate passed) OR documented divergence (gate failed → reverted).

**Expectation:** any-overlap vs reciprocal-50% are definitionally different; the catalog run will likely change which copies collapse → gate FAIL on `cat.families.tsv`/`cat.copies.tsv`. Expected; document and defer.

- [ ] **Step 1: Attempt the routing**

Replace the two `distinct_locus_reps(...)` calls in denovo_pipeline.rs with the reciprocal-50% predicate (`distinct_loci` logic from driver.rs:212 — adapt its signature to the `DenovoTranscript` copies these call sites pass). Keep `distinct_locus_reps` defined.

- [ ] **Step 2: Build + test**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -5
cargo test --lib 2>&1 | tail -5
```
Expected: builds; some `distinct_locus_reps_*` tests may now describe the old predicate — note them, do not delete.

- [ ] **Step 3: Gate check**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/mechanism/byte_identity_gate.sh check
```
- **If GATE PASS:** the predicates agreed on the corpus — commit as landed (as in Task 5 Step 3), and update the two tests to the retained predicate.
- **If GATE FAIL:** Step 4.

- [ ] **Step 4 (if gate failed): revert, document, defer**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git checkout -- src/rustle/vg_family/denovo_pipeline.rs
cargo build --release 2>&1 | tail -3
bash bench/mechanism/byte_identity_gate.sh check   # PASS again
```
Append to `consolidation_divergences.md`:
```
| B2.3 >=2-loci | distinct_loci (reciprocal-50%) | cat.families/copies differ: any-overlap (distinct_locus_reps, live) collapses loci the reciprocal-50% rule keeps separate (or vice-versa) — definitionally different predicates. | deliverable C |
```
Commit the doc:
```bash
git add bench/mechanism/consolidation_divergences.md
git commit -m "docs(consolidation): B2.3 >=2-loci predicates DIVERGE (any-overlap vs recip-50%) -> defer to C

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: B2.4 — attempt refiner consolidation (CNM vs Louvain)

**Files:**
- Reference: `src/rustle/vg_family/family_definition.rs:437` `refine_families` (CNM), `src/rustle/vg_family/family_split.rs:477` `gamma_quasi_clique_partition` (Louvain)
- Modify: whichever caller sits on the live path (trace it in Step 1)

**Interfaces:**
- Produces: EITHER landed (gate passed) OR documented divergence (gate failed → reverted).

**Expectation:** two different community-detection algorithms; likely divergent on the catalog run. Expected; document and defer.

- [ ] **Step 1: Trace which refiner is live and who calls the other**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
grep -rn "refine_families\b\|gamma_quasi_clique_partition" src/ | grep -viE "test|fn refine_families\b|fn gamma_quasi" | head
```
Expected: identify the live caller (per ONE_METHOD, `gamma_quasi_clique_partition`/Louvain is on the live path; `refine_families`/CNM is parity-tested). The consolidation target is the parity-tested one — route the live caller to it.

- [ ] **Step 2: Attempt routing the live caller to the parity-tested refiner**

Replace the live refiner call with the canonical (parity-tested) one, matching signatures. Keep both fns defined.

- [ ] **Step 3: Build + test + gate**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -5
cargo test --lib 2>&1 | tail -5
bash bench/mechanism/byte_identity_gate.sh check
```
- **If GATE PASS:** commit as landed (retain the canonical refiner; note the other is now unreferenced-but-kept for compile-pinning).
- **If GATE FAIL:** Step 4.

- [ ] **Step 4 (if gate failed): revert, document, defer**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git checkout -- src/rustle/vg_family/   # revert the routing edit
cargo build --release 2>&1 | tail -3
bash bench/mechanism/byte_identity_gate.sh check   # PASS again
```
Append to `consolidation_divergences.md`:
```
| B2.4 refiner | parity-tested (CNM) | cat.families differ: CNM and Louvain partition the homology graph differently on real families — the two are not algorithmically equivalent. | deliverable C |
```
Commit:
```bash
git add bench/mechanism/consolidation_divergences.md
git commit -m "docs(consolidation): B2.4 refiners DIVERGE (CNM vs Louvain partitions) -> defer to C

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: Update the Consistency table + final gate

**Files:**
- Modify: `docs/ONE_METHOD.md` (the "Consistency status" table)

- [ ] **Step 1: Reflect actual outcomes in the table**

For each B2 item, update its row in `docs/ONE_METHOD.md`'s Consistency table to state the real outcome: **consolidated (byte-identical, gate passed)** or **diverges → deferred to deliverable C (see consolidation_divergences.md)**. Correct the table's blanket "none of these change the numbers" line to: "the ones that consolidated are byte-identical (gated); the ones that diverged are documented — they were never equivalent, and the correct choice is deliverable C."

- [ ] **Step 2: Final full gate + cargo test**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/mechanism/byte_identity_gate.sh check
cargo test --lib 2>&1 | tail -5
grep -rn '"/home/juanfra/winloci_scratch' src/ | grep -viE 'test|///|//!' || echo "no silent scratch defaults"
```
Expected: `GATE PASS` (B's invariant: baseline reproduced at B's end); tests green; no silent scratch defaults.

- [ ] **Step 3: Commit**

```bash
git add docs/ONE_METHOD.md
git commit -m "docs(consolidation): update Consistency table with B outcomes (consolidated vs deferred-to-C)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review Notes (for the planner)

- **Spec coverage:** §2 gate → Task 1. §3 B1a paths → Task 2; B1b flags → Task 3. §3 B2 items 1–4 → Tasks 4–7. §4 escalation → the Step-4 branch in Tasks 5–7 + `consolidation_divergences.md`. §6 success criteria → Task 8. All covered.
- **The invariant is enforced mechanically:** every task ends with `byte_identity_gate.sh check` PASS — either because the change was byte-identical, or because a diverging change was reverted. No task can leave the tree with a changed number.
- **Divergence is a first-class outcome:** Tasks 5–7 are structured so "gate fails → document → defer" is a *complete* task, not a failure. The plan explicitly expects EM / ≥2-loci / refiner to diverge (different error model, different predicate, different algorithm) and lands them as findings.
- **Placeholder scan:** the routing edits in Tasks 5–7 describe the signature mapping rather than showing final code, because the exact adaptation depends on the live call-site's local variables — but each names the exact functions, files, line numbers, and the decision (gate PASS/FAIL branch) concretely. No thresholds or types are left unnamed.
- **No new numbers introduced:** B adds no constant; it only removes duplication or exposes existing constants as flags with identical defaults.
