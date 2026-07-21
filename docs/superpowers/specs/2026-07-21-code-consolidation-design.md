# Deliverable B — Byte-Identical Code Consolidation — Design

**Date:** 2026-07-21
**Status:** approved design, pre-implementation
**Predecessor:** Deliverable A (mechanism transparency, `2026-07-20-mechanism-transparency-design.md`).
A's artifact (`bench/rustle_mechanism.html`) discloses that several mechanisms have duplicate
implementations coexisting in the code and points here for the fix. B removes the duplication
**without changing any scientific number.**

**Motivation.** The advisor's "you jump between approaches" reading is fed by the code literally
carrying two refiners, two "≥2 loci" predicates, and three EM paths (self-documented in
`docs/ONE_METHOD.md`'s "Consistency status" table). A made the single method *legible*; B makes
it *singular* in the code, so the artifact's "these are stages of one pipeline" claim is backed
by there being one implementation to describe.

---

## 1. The hard invariant

**B never changes a number.** Every consolidation that lands leaves the real-data outputs
**md5-identical** to the pre-B baseline. Safety is **proven per item** by an empirical
byte-identity gate — never assumed from the docstring claim that "these don't change the
numbers" (that claim is a hypothesis per item; see §4's escalation rule for when it fails).

Two things are explicitly **out of scope** for B because they *would* change numbers:
- The five inconsistent per-base error rates (`error_rate=0.003`, `QUANT_ERROR=0.01`,
  `MOSAIC_EPS=0.01`, mosaic `eps_floor=0.005`/`eps_cap=0.05`, `RUSTLE_CONFLICT_EPS=0.001`) →
  **deliverable C**.
- The assembler carve (41 modules, compile-pinned by ~12 tendril symbols) → a separate effort
  (`RETIREMENT_AND_MIGRATION.md`).

---

## 2. The byte-identity gate (built first, before any consolidation)

A single harness that is the test for every change in B.

**`bench/mechanism/byte_identity_gate.sh`** — runs a fixed corpus on the real gorilla data
(`/home/juanfra/winloci_scratch/GGO_mm.bam` + `GGO.fasta`, foreground per the crash rule) and
writes md5s of the output TSVs to `bench/mechanism/byte_identity_baseline.txt`:

| Run | Command (region-restricted, foreground) | md5'd outputs |
|---|---|---|
| GSTM assignment | `copy_assign --region NC_073224.2:129160000-129230000 --homology-primary --min-copies 2 --dump-psv` | `.assignments.tsv`, `.families.tsv`, `.quant.tsv`, `.psv_cols.tsv` |
| Catalog (small) | `gw_family_catalog` on a fixed 2–3 region list | `.families.tsv`, `.copies.tsv` |
| ASJ | `asj --region <ABCC4>` | `.asj.tsv` |

- The gate has two modes: `freeze` (write baseline) and `check` (re-run, diff against baseline,
  exit non-zero on any mismatch, naming the drifted file).
- **Baseline frozen once at B's start** (current HEAD, commit recorded in the baseline file).
- Every consolidation runs `check`; **PASS iff every md5 equals baseline.**
- The EM consolidation additionally freezes `--em` / `--vg-realign` / default `quant.tsv` md5s,
  since EM is the leg most likely to diverge.

This reuses the GSTM-md5 discipline proven in the vg-o2-substrate work (where
`.assignments.tsv` md5 `8c0507aa` gated a representational refactor).

---

## 3. The two phases

### B1 — config hygiene (byte-identical by construction)

**B1a — fix the 6 hardcoded `winloci_scratch` paths.** They are compile-time string defaults
that silently point at one machine:
- `asj_strand_bias.rs:58` `DEFAULT_BAM`
- `asj_verify.rs:75` `DEFAULT_BAM`
- `driver.rs:103` `SCRATCH`
- `o2_materialize.rs:59` `FAM_TSV_DEFAULT`, `:255` `BAM_DEFAULT`, `:421` `GENOME` default

Action: remove the hardcoded default; make the value a **required CLI arg / explicit parameter**
(or error-if-unset), never a silent fallback to `/home/juanfra/...`. Byte-identical because the
thesis runs always pass explicit paths — the gate confirms it. A machine without those paths now
gets a clear error instead of a wrong silent default (a latent-bug fix, disclosed).

**B1b — expose result-guards + buried decision constants as CLI flags, defaults unchanged.**
- `POA_CAP` (copy_assign_pipeline.rs:524), `READ_CAP` (o2_materialize.rs:250) → `--poa-cap`,
  `--read-cap` (or env), default = current value → byte-identical by construction.
- Any decision-tier constant currently reachable only by editing source and not already a flag
  (audit against A's `heuristics.tsv` decision rows vs the existing CLI surface).
- **Not** the 137 inert constants — exposing them all is YAGNI and clutters the CLI.

### B2 — code consolidations (each md5-gated, one at a time)

Per `docs/ONE_METHOD.md`'s Consistency table, in this order (lowest-divergence-risk first):

| # | Consolidation | Canonical target | Divergence risk |
|---|---|---|---|
| 1 | `family_define` binary → `--legacy` parity fixture | `gw_family_catalog` is sole catalog | low (already parallel) |
| 2 | Three EM paths → one | `em_copy_assign::em_assign_family` | medium (verify `quant.tsv`) |
| 3 | Two "≥2 distinct loci" predicates → one | reciprocal-50% (`distinct_loci`) | **high** (any-overlap is looser) |
| 4 | Two γ-quasi-clique refiners → one | the parity-tested one | **high** (CNM vs Louvain) |

Each item follows the **protocol**:
1. Confirm the gate baseline is frozen.
2. Route all callers of the redundant impl to the canonical one (smallest diff; leave the
   redundant fn in place but unreferenced, or feature-gate it, if deletion pins compilation).
3. `cargo build` + `cargo test` green.
4. Run the byte-identity gate `check`.
5. **If all md5 identical →** commit as byte-identical-safe (message states the gate passed).
   **If any md5 differs →** the escalation rule (§4).

Items 3 and 4 are flagged high-risk precisely because the two impls are *definitionally* or
*algorithmically* different; the gate is expected to catch a real divergence there, and that is
a successful outcome of the protocol, not a failure.

---

## 4. The escalation rule (what happens when the gate fails)

When a consolidation changes any md5:
1. **Do not merge it.** B's invariant holds — outputs stay byte-identical.
2. **Keep both implementations.** Revert the routing; the redundant impl stays.
3. **Document the divergence as a finding:** which families/regions differ, by how much
   (family count, copy count, or the specific rows), written to
   `bench/mechanism/consolidation_divergences.md`.
4. **Surface it.** A divergence means the two implementations were *never equivalent* — the
   docstring's "doesn't change the numbers" was wrong for that item. That is defensible data for
   the advisor (and a candidate for deliverable C, where the *correct* impl is chosen with a
   before/after).

So B produces two kinds of result: consolidations that landed (duplication removed, numbers
identical) and divergences that were caught and documented (duplication retained, finding
recorded). Both are honest outcomes.

---

## 5. Components & boundaries

| Unit | Purpose | Depends on |
|---|---|---|
| `bench/mechanism/byte_identity_gate.sh` | freeze/check real-data md5 baseline | built binaries, GGO data, crash rule |
| `bench/mechanism/byte_identity_baseline.txt` | the frozen md5s + base commit | the gate |
| `bench/mechanism/consolidation_divergences.md` | log of caught divergences (may be empty) | escalation rule |
| the routing edits | point callers at one impl each | the gate as their test |

Each consolidation is independently reviewable: a reviewer can reject item 3 while approving
item 1, and each carries its own gate-pass evidence.

---

## 6. Success criteria

- The byte-identity gate exists, is frozen at B's start, and every landed commit passes `check`.
- B1a: no `/home/juanfra/winloci_scratch` string default remains as a silent fallback in any
  shipped binary path (grep-clean); the gate still passes.
- B1b: `POA_CAP`/`READ_CAP` (and any buried decision constant) are CLI-settable; defaults
  unchanged; gate passes.
- B2: each of the 4 consolidations either (a) landed with a gate-pass, or (b) is documented in
  `consolidation_divergences.md` with the specific divergence and left un-merged.
- `cargo test` green throughout.
- `docs/ONE_METHOD.md`'s Consistency table updated to reflect what consolidated vs what diverged.
- **No scientific number changed** (the invariant), verifiable by the gate baseline being
  reproduced at B's end.

---

## 7. Risks

- **A consolidation deletes a fn that pins compilation** (the "12 tendril symbols" problem).
  Mitigation: route callers without deleting; feature-gate or `#[allow(dead_code)]` the
  redundant impl rather than removing it, unless removal is clean.
- **The gate corpus is too small to catch a rare divergence.** Mitigation: the corpus targets
  the exact mechanisms being consolidated (GSTM exercises assignment+EM; the catalog run
  exercises refiner+≥2-loci); expand the region list if a consolidation touches an untested path.
- **Items 3/4 diverge and B lands fewer consolidations than hoped.** That is the honest, correct
  outcome — a divergence is a finding, not a failure, and the invariant is preserved either way.
