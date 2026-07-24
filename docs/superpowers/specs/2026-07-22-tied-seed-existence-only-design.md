# Tied-seed copies as existence-only appends — design

**Date:** 2026-07-22
**Status:** approved design, pre-implementation
**Type:** correctness fix / restructure of `detect_and_assign` (opt-in path only)

## Motivation

`--tied-seed` (opt-in) currently extends the de-novo `skeletons` with tied-seed skeletons before
`assemble_gate` + `collapse_loci_span_aware`, so the tied-seed transcripts enter the same `reps` set as
the primary transcripts. Those tied reps then flow through the ENTIRE downstream structuring pipeline —
`build_read_placements` → `conflict_edges` → `conflict_families` → `refine_families_exon_sum` →
`assign_family_detailed`. Because tied-seed copies are seeded from K=0 (MAPQ-0) reads, in the conflict
graph they tie with everything, adding spurious edges that OVER-MERGE families and wreck assignment.

Measured on the chr1 amylase family (ID_131, `copy_assign --region chr1:103415755-103845829`):

| stage             | baseline (no tied-seed) | `--tied-seed` |
|-------------------|-------------------------|---------------|
| conflict edges    | 83                      | **107**       |
| families (refine) | 3                       | 2             |
| **copies (quant)**| **21**                  | **6**         |
| read assignments  | 94,195                  | **27,187**    |

The tied-seed run DROPS the real annotated amylase genes AMY1A/1B/1C/AMY2A (21→6 copies) and craters
assignment. The genome-wide sensitivity sweep confirmed the net effect is NEGATIVE: over the 9 co-located
Soto families with missed members, `--tied-seed` recovered 2 members (TRIM64B, TRIM64) but destroyed 4
(the amylase genes), net −2. The same joint-collapse mechanism produced the os1 rep-shift artifact.

An earlier attempt to fix this by collapsing primaries separately and appending tied reps to `reps`
additively FAILED (amylase still 6): the tied reps still enter `reps` and poison conflict/refine/assign
wherever they sit. The tied reps must be kept OUT of the structuring pipeline entirely.

## Design

Tied-seed copies become **existence-only appends**, mirroring the collapse-gate block that already
exists at the tail of `detect_and_assign` (which appends `gated_family(rep, …)` for collapsed reps not in
any assembled family).

### Flow (when `cfg.tied_seed`)

1. **Primary pipeline unchanged.** `reps` = primary skeletons only (do NOT `skeletons.extend(tied)`).
   `build_read_placements` / conflict / refine / assignment / collapse-gate all run exactly as in a
   baseline run, producing `out: Vec<FamilyAssignment>` byte-identical to today's no-`--tied-seed` output.
2. **Tied reps computed separately.** `tied_seed_skeletons(rescue_extra, &skeletons, cfg.pass1_min_reads)`
   → `assemble_gate(&tied_sk, genome, &cfg.gate)` → (readthrough/mischain filter, same as primary) →
   `collapse_loci_span_aware` → `tied_reps: Vec<DenovoTranscript>`. These NEVER enter `reps`.
3. **Existence-only append.** After `out` is fully built, gather every emitted copy span from
   `out[*].copy_spans`. For each tied rep whose `(chrom, start, end)` overlaps NO emitted copy span,
   append an existence-only family built like `gated_family`: every overlapping non-supplementary read
   marked `AssignStatus::Tied` (unassignable, `posterior = [1.0]`), `n_copies = 1`, `collapsed_copies`
   recorded, `family_id = "TSFAM{out.len()}"`, `copy_tids/copy_spans/copy_introns` = the rep's.

### Grouping decision

Each admitted tied rep is its own **singleton existence-only family** (`n_copies = 1`). Rationale: the
Soto member-DETECTION metric only needs a copy overlapping the member; attaching tied copies to a
homologous primary family would re-introduce the homology + merge machinery this fix exists to bypass,
for no detection benefit. (A tied copy is conceptually a member of its tie partner's family, but that is
an O2/assignment relationship we deliberately abstain on.)

### Overlap criterion

A tied rep is DROPPED if its span intersects any emitted copy span on the same chrom
(`r.start < t.end && t.start < r.end`) — the same span-intersect test used elsewhere. This keeps
tied-seed strictly additive: it can only add copies at loci no primary/collapsed copy already occupies.

### New helper (unit-testable, pure)

```
fn tied_existence_families(
    tied_reps: &[DenovoTranscript],
    emitted_spans: &[(String, u64, u64)],   // all copy spans already in `out`
    bam_reads: &[BamRead],
    next_family_index: usize,
) -> Vec<FamilyAssignment>
```
Returns one existence-only `FamilyAssignment` per non-overlapping tied rep, family ids
`TSFAM{next_family_index + k}`, every overlapping read `AssignStatus::Tied`. Deterministic (input order).

### Marking / transparency

The internal `TSFAM` family-id is for readability inside `detect_and_assign` only — the `copy_assign`
emitter RENUMBERS every returned family to `CAFAM{gfam}` (copy_assign.rs:~1034), so the `TSFAM` prefix
does NOT reach `family.tsv`/`quant.tsv`. In the emitted output a tied-seed existence copy is identifiable
by its SIGNATURE, not a label: `n_copies == 1`, every read `AssignStatus::Tied`, and `abundance == 0.0`
(no assigned reads). The stderr log carries the count:
`[detect_and_assign] tied-seed existence-only: N tied reps -> M appended (non-overlapping)`.

DEFERRED (follow-up, not this change): if analysts need to filter tied-seed existence copies directly, add
an explicit provenance marker column (or propagate the `Skeleton.tied_seeded` flag) to the quant row —
today the honest K=0-existence-vs-K≥1-assigned split must be reconstructed from the signature above.

## Byte-identity guarantee

`--tied-seed` OFF: the `if cfg.tied_seed` block is skipped entirely; `skeletons` is never extended and no
existence families are appended → byte-for-byte identical to today. The md5 byte-identity gate must pass
for the OFF path.

## Acceptance criteria (TDD targets)

1. **Off = identical.** `DenovoConfig::default().tied_seed == false`; a detect_and_assign run with
   `tied_seed = false` produces the same families/copies/assignments as before this change.
2. **`tied_existence_families` unit test.** Given emitted spans and tied reps, only non-overlapping tied
   reps become families; each is `n_copies = 1`, `AssignStatus::Tied`, `TSFAM*` id; overlapping tied reps
   are dropped.
3. **Amylase preserved (integration, verified via `copy_assign` rerun).** `--tied-seed` on
   chr1:103415755-103845829 → 21 copies, AMY1A/1B/1C/AMY2A all overlapped, primary assignment count ==
   baseline (94,195).
4. **TRIM64 recovered (integration).** `--tied-seed` on chr11:89646664-90033077 → the 11 baseline primary
   copies unchanged + TRIM64B (chr11:89,785,154) and TRIM64 (chr11:89,893,563) appended as `TSFAM`
   existence families.
5. **Appended copies are unassignable.** Every read in a `TSFAM` family has `AssignStatus::Tied`.

## Out of scope

- Genome-wide promotion (`gw_family_catalog` still has no tied-seed) — a separate follow-up once this is
  proven safe per-region.
- Any change to the primary conflict/refine/assignment logic.
- Attaching tied copies to homologous families (deferred; singleton existence families for now).
- The `de`-tie gate and `RUSTLE_TIED_SEED_DEBUG` diagnostic (already in the tree, unaffected).
