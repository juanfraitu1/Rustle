# Tie-break invariance certificate — Design

**Date:** 2026-07-10. **Scope:** narrow, additive, report-only. **Motivation (advisor).** The primary/secondary
alignment label is arbitrary at MAPQ = 0, so a locus defined by "≥ 1 primary read" is only sound if its existence
is invariant under that arbitrary choice. The invariance experiment
(`bench/PRIMARY_SECONDARY_INVARIANCE.md`) measured this per family; this design makes the tool **emit the
certificate itself**, per copy, so every recovered copy carries whether its existence survives relabeling.

## What it computes

For each copy `ci` of a family:

- **anchored_reads** = number of assigned reads whose primary alignment maps *uniquely* (`mapq > 0`) and whose
  `best_copy == ci`. minimap2 sets `mapq = 0` exactly when the primary is tied/arbitrary, so a `mapq > 0` read is
  NOT a candidate for relabeling — its support for `ci` is fixed under every tie-break. This reproduces the
  experiment's `samtools view -c -F 2308 -q 1 <copy_span>` bound, computed in-binary from the reads already loaded.
- **tie_invariant** = `anchored_reads ≥ GATE_MIN_READS` (= 3, the same gate a locus must clear to exist,
  `denovo_assemble::GATE_MIN_READS`). TRUE ⟹ the copy is supported by ≥ 3 unique mappers ⟹ it exists under EVERY
  primary/secondary relabeling (adversarially proven). FALSE ⟹ fewer than 3 unique mappers ⟹ its existence relies
  on the arbitrary primary assignment of tied reads and/or on junction structure — flagged for scrutiny.

**Conservative by design.** `tie_invariant` is the *unique-mapper* bound only. A junction-defined copy (DAZ2: 1
unique mapper, but recovered via its 31-vs-16 intron structure and shown invariant by the end-to-end flip) reads
FALSE here. FALSE therefore means "not guaranteed by unique support alone," NOT "spurious." This is stated in the
column doc and the bench note so the flag is not over-read.

## Output

`quant.tsv` gains two trailing columns:

```
family_id  copy_index  copy_tid  copy_chrom  copy_start  copy_end  abundance  ci95_halfwidth  n_reads_hard  anchored_reads  tie_invariant
```

Recovery is **unchanged** — no copy is added, dropped, or reassigned; only two descriptive columns appear. A
one-line stderr summary reports `N/M copies tie-break-invariant`.

## Non-goals

- **No recovery change.** This is the report-only option; it does not gate, abstain, or re-assign (those were the
  rejected "honest abstain" / "invariant support gate" alternatives).
- **No junction-invariance signal in v1.** The certificate is the unique-mapper bound; junction-based invariance
  (what actually rescues DAZ2) is demonstrated by the flip experiment, not computed here. Noted as a future column.

## Files

- **Modify** `src/bin/copy_assign.rs`: add `anchored: usize` + `tie_invariant: bool` to `QuantRow`; a pure
  `anchored_support(best_copies, mapqs, ci) -> usize` helper (unit-tested); populate at the existing per-copy
  loop (`fa.assignments` + `bam_reads[ri].mapq` in scope); extend the header + writeln; add the stderr summary.
- **Modify** `bench/PRIMARY_SECONDARY_INVARIANCE.md`: note the certificate is now emitted, with the
  conservative-bound caveat.

## Testing (TDD)

- `anchored_support`: `best=[0,0,1,0]`, `mapq=[60,0,60,60]`, `ci=0` ⟹ 2 (the mapq=0 read excluded); `ci=1` ⟹ 1.
- `anchored_support` with all `mapq=0` ⟹ 0 for every copy (the TSPY case).
- `tie_invariant` boundary: anchored 2 ⟹ false, 3 ⟹ true (`>= GATE_MIN_READS`).

## Reproduce / expected

```
copy_assign ... --region NC_073248.2:42778133-42950552 --out kf_DAZ
# quant.tsv: DAZ1 anchored~178 tie_invariant=TRUE ; DAZ2 anchored~1 tie_invariant=FALSE (junction-defined)
copy_assign ... --region NC_073248.2:34731504-34847734 --out kf_TSPY
# quant.tsv: all copies anchored=0 tie_invariant=FALSE (the tie-break-dependent array)
copy_assign ... --region NC_073224.2:129160000-129240000 --out kf_GSTM
# quant.tsv: all copies tie_invariant=TRUE (unique mappers)
```

Related: `bench/PRIMARY_SECONDARY_INVARIANCE.md`, `project_primary_secondary_invariance`,
`project_k0_frontier_unresolvable` (TSPY = the frontier), `denovo_assemble::GATE_MIN_READS`.
