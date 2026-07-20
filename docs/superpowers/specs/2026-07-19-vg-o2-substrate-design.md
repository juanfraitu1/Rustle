# The Family Variation Graph as the O2 Assignment Substrate — Design Spec

**Status:** approved (brainstorming), ready for implementation plan
**Date:** 2026-07-19
**Builds on:** the O2 copy-assignment decision (`copy_assign.rs`: `read_copy_evidence` path-log-likelihood + `assign_read_editing` significance gate) and the copy-graph emitter (`copy_graph.rs`: read `W`-lines, copy `P`-lines).
**Objective:** make the per-family variation graph the **actual substrate the O2 decision runs on** — an **auditable ad-hoc reference** (copies = paths, built from the data, not the linear genome) through which each read is **threaded as a walk whose identity and assignment stay traceable** — while keeping the decision **bit-identical** and the significance certificate untouched.

## Goal

Convert the O2 per-read decision from an implicit `Vec<Option<u8>>` matrix into an explicit `BubbleGraph` (PSV bubbles + copy-paths) that the decision constructs and threads reads through, and emit that decision graph with per-read provenance (walk + assigned copy + significance certificate). `.assignments.tsv` stays **byte-identical**; the change is representational + additive.

## Motivation

The advisor's standing concern is **reference bias**: a single linear reference silently drops or mis-places reads that do not match it, and collapses read-level evidence into reference-anchored counts. The graph-fraction audit (`wf_f5e1a2f6`) found the O2 **decision** to be the pipeline's biggest "graph" honesty gap: it is mathematically a per-copy path-log-likelihood, but it runs on a matrix and argmaxes over `K` scalar path-scores — no graph object is constructed or traversed, and the materialized VG (`assemble_vg`, `copy_graph`) is only *emitted*, never *consumed* by the decision. This change closes that gap honestly: the **ad-hoc reference** is the family's own copies as graph paths (no linear allele privileged), the decision **threads** each read through it, and the graph + read-walks + per-read certificate are emitted as the audit trail — so "assignment = threading a read through the variation graph" and "which reads make each copy, auditably" become literally true, not documented frames. The decision numbers do not change (the arithmetic is identical), so no O2 validation re-opens.

## Design

### §1 — The ad-hoc reference: an explicit `BubbleGraph`

A new type (in `copy_assign.rs` or a sibling `bubble_graph.rs`) built **once per family** from `&[CopyProfile]`:
```rust
pub struct BubbleGraph { pub bubbles: Vec<Bubble>, pub n_copies: usize }
pub struct Bubble {
    pub col: usize,                    // PSV column index (= bubble id); SAME order as the matrix
    pub copy_allele: Vec<Option<u8>>,  // the allele-node each copy PATH visits here (index = copy); None = gap
    pub decisive: bool,                // copies carry >=2 distinct non-None alleles here (read-independent)
}
```
`copy_allele[ci] = copies[ci].alleles[col]`; `decisive` = the distinct non-`None` alleles number ≥ 2. Each **copy is a path** through the bubbles (its allele-node at each), plus its copy-specific junctions (`CopyProfile.junctions`) as path edges. This IS the ad-hoc reference — the family's copies, built from the data, per family, auditable. Constructed via `BubbleGraph::from_copies(&[CopyProfile])`.

### §2 — The decision threads reads through it (byte-identical)

Refactor `read_copy_evidence` (`copy_assign.rs:267`) to take the caller-prebuilt `&BubbleGraph` (for the log-likelihood walk) **plus** the `&[CopyProfile]` slice it already receives (the junction term and the `copy_pair_significance` certificate still read the copies directly, so they stay bit-identical). It then **threads the read as a walk**: iterate the read's spanned bubbles (`read.psv_obs[b.col].is_some()`), and for each, accumulate the per-copy path term exactly as today —
```
lp_match = (1-e).ln(); lp_mis = (e/3).ln();     // e from psv_qual[col] or p.error_rate
if b.decisive { n_decisive += 1 }
for ci in 0..n { if let Some(a) = b.copy_allele[ci] { logl[ci] += if obs==a {lp_match} else {lp_mis} } }
```
The junction term, the `best = argmax(logl)` (earliest on ties), and `min_p` (the IsoCon identifiability bound via `copy_pair_significance`) are computed exactly as now. **The f64 terms are summed in the same order over the same columns**, so `logl`, `min_p`, `n_decisive` are bit-identical. `decisive` is precomputed once per family (read-independent) instead of re-derived per read — same values, and a mild perf win (the current code re-scans all copies per column per read).

### §3 — The significance gate is untouched

`assign_read_editing` (`copy_assign.rs:367`) consumes the `ReadEvidence { logl, min_p, n_decisive }` unchanged: the same argmax + runner-up + margin + Bonferroni `alpha/(n-1)` certificate → `Assigned` / `Ambiguous` / `Tied`, and the same softmax posterior. The EM path (`read_copy_evidence` is shared) is likewise unchanged. **No threshold, no p-value, no `1/k`, no gate logic is modified** — the certificate remains the anti-arbitrary-threshold virtue.

### §4 — Per-read provenance + audit emit (the auditable trail)

The decision retains each read's **walk** (`ReadFeatures.psv_obs` over the bubbles) + **assigned copy** + **certificate** (`p_read`, `min_p`, `status`, `n_decisive`) — never collapsed to a bare per-copy count. Emit the **decision graph** as an auditable GFA (`<out>.vg_audit.gfa`), extending the existing `copy_graph.rs` emitter:
- copy `P`-lines = the ad-hoc reference (each copy's path through the bubble allele-nodes);
- read `W`-lines (already emitted via `ReadWalk`, `copy_graph.rs:76/225`) **tagged** with the per-read certificate: `CP:Z:<assigned_copy> PV:f:<p_read> MP:f:<min_p> ST:Z:<Assigned|Ambiguous|Tied>`.
So a reader can inspect, on the exact graph the decision used, which reads make each copy and with what confidence. The emit is **opt-in** (a flag, e.g. `--vg-audit`; folded into `--phase` if that is cleaner) so the default `.assignments.tsv` run is byte-identical and un-slowed.

### §5 — Byte-identity is the safety property

`.assignments.tsv` (and the EM quant) must be **bit-identical** to the pre-change output on every family. This is provable by construction (identical arithmetic/order) and checked by (a) the existing `read_copy_evidence_matches_assignment_internals` test and (b) a parity test that asserts the threaded `ReadEvidence` (`logl`/`min_p`/`n_decisive`) equals the **golden pre-change values** — captured from the current implementation before the refactor and frozen as fixtures — on unit cases, plus a **byte-identical `.assignments.tsv`** diff on ≥1 real family. Because the decision is unchanged, no downstream O2 validation (SUN determinism, significance-gate rates, etc.) re-opens.

### §6 — Out of scope

- Changing any assignment, the significance gate, the softmax posterior, or the EM.
- O1 family definition (untouched).
- **Reference-bias RECOVERY** — re-threading reads the *linear* reference lost/mis-placed against the ad-hoc reference (the deferred second stage the user flagged as "next"). This spec makes the decision run on, and be auditable via, the graph; it does **not** re-map reads.
- Unifying with the Python-parity `assemble_vg` emission (its own path stays; the audit emit comes from the decision's `BubbleGraph`).

### §7 — Testing

- **Parity (the safety gate):** `BubbleGraph::from_copies` correctness (bubble order, `copy_allele`, `decisive`); the threaded `read_copy_evidence` returns a `ReadEvidence` bit-identical to the **golden captured values** on unit fixtures (SNP-bubble, gap, junction, editing-column); and a data-gated real family (`copy_assign` foreground/serial/`winloci_scratch`, crash rule) with `.assignments.tsv` unchanged byte-for-byte.
- **Audit emit:** the `<out>.vg_audit.gfa` is valid (every `P`/`W` step backed by an `L`-line, as `copy_graph`'s existing invariant test checks), copy paths match the ad-hoc reference, and each read `W`-line carries the correct `CP/PV/MP/ST` tags matching its `.assignments.tsv` row.

## Global constraints

- **Bit-identical `.assignments.tsv`** (and EM quant): the decision arithmetic and column order are preserved exactly; the graph is a re-expression, not a new algorithm.
- The significance certificate / gate is **not modified**.
- The `BubbleGraph` is built once per family; `read_copy_evidence` and the EM share it.
- The audit emit is **opt-in** and additive; default runs are byte-identical and not slowed.
- No O1, no mapping/recovery, no `assemble_vg` unification.
- Validation runs of `copy_assign` are **foreground + serial + `winloci_scratch`** (WSL2 crash rule); no `copy_assign` background/nohup/pkill.
