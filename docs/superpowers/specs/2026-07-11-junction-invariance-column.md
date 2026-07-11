# Junction-invariance column — Design

**Date:** 2026-07-11. **Scope:** narrow, additive. Extends the tie-break invariance certificate
(`2026-07-10-tie-invariance-certificate.md`) with the missing mechanism: a copy identifiable by its
**copy-specific junctions** is invariant to the arbitrary primary/secondary label even when it has no unique
mappers. This is what actually rescues DAZ2 (1 unique mapper → `tie_invariant=FALSE`, yet recovered via its
31-vs-16 intron structure and shown invariant by the end-to-end relabel experiment).

## What it computes

A read is **junction-only** support for its assigned copy when a copy-specific junction made it decisive but
PSVs alone could not — the existing `junction_only` predicate (`r.combined.n_decisive >= 1 && r.psv.n_decisive
== 0`). Such a read carries a junction distinctive to that copy, so it identifies the copy by splice structure
**regardless of which locus holds its primary alignment** — a relabeling-invariant signal.

- **copy_junction_support[ci]** = number of junction-only reads whose `best_copy == ci` (new per-copy tally,
  bucketed at fa-build time where `r.psv` and `r.combined` are both in scope).
- **junction_invariant** (new `quant.tsv` column) = `copy_junction_support[ci] >= GATE_MIN_READS` (3).

**Correctly discriminating.** A single-exon readthrough (RFPL's artifacts) carries no junctions ⟹ 0 junction
support ⟹ FALSE — the column does not bless artifacts. Exonically-identical copies (TSPY) share all junctions ⟹
no copy-specific junction is decisive ⟹ FALSE — correctly the K=0 frontier. A junction-distinguished copy
(DAZ2) ⟹ TRUE.

## Output

`quant.tsv` gains one trailing column after `tie_invariant`:

```
... n_reads_hard  anchored_reads  tie_invariant  junction_invariant
```

The stderr summary is updated to the **overall** bottom line: a copy is invariant if `tie_invariant OR
junction_invariant` (unique mappers OR distinctive junctions), so DAZ2 now counts as invariant.

Recovery is **unchanged** — additive columns only.

## Files

- **Modify** `src/rustle/vg_family/denovo_pipeline.rs`: add `pub copy_junction_support: Vec<usize>` to
  `FamilyAssignment`; init `Vec::new()` in `empty()` and `vec![0; all_copies.len()]` at the main build; a pure
  `add_junction_support(support, best_copy, combined_decisive, psv_decisive)` helper called in the per-read
  fa-build loop (mirrors the existing `junction_only` line); the compiler flags any test construction sites
  (add `Vec::new()`).
- **Modify** `src/bin/copy_assign.rs`: add `junction_invariant: bool` to `QuantRow`, populate from
  `fa.copy_junction_support`, extend header + writeln, update the stderr invariance summary to the
  `tie_invariant OR junction_invariant` count.

## Testing (TDD)

- `add_junction_support`: combined decisive + PSV NOT decisive, `best_copy=1`, `[0,0]` ⟹ `[0,1]`; combined +
  PSV both decisive ⟹ unchanged (PSV resolved it); combined NOT decisive ⟹ unchanged; `best_copy` out of range
  ⟹ unchanged (guarded).
- Existing `junction_resolves_read_when_psvs_cannot` pins the underlying per-read predicate.

## Reproduce / expected

```
copy_assign ... --region NC_073248.2:42778133-42950552 --out kf_DAZ
# quant.tsv: DAZ2 tie_invariant=false BUT junction_invariant=TRUE (its distinctive introns)
copy_assign ... --region NC_073248.2:34731504-34847734 --out kf_TSPY
# quant.tsv: junction_invariant=FALSE for all copies (exonically identical -> no copy-specific junction)
copy_assign ... --region NC_086018.1:30200000-30390000 --out kf_RFPL
# quant.tsv: the readthrough copies junction_invariant=FALSE (single-exon, no junctions)
```

Related: `docs/superpowers/specs/2026-07-10-tie-invariance-certificate.md`,
`bench/PRIMARY_SECONDARY_INVARIANCE.md`, `project_primary_secondary_invariance`,
`project_daz2_locus_support` (DAZ2 = junction-defined).
