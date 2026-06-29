# Multimapper rescue report — the "missing side" as a usable output

The genome-wide analysis quantified the treasure StringTie discards: ~424 annotated transcripts
recoverable *only* from multimappers, plus the disambiguation evidence that separates them from
~456 phantom copies. This feature turns that into a **first-class output**: a TSV, written next
to the GTF, that lists per VG-family transcript what a primary-only assembler sees versus the
decisive evidence recoverable from multimappers — with a `rescued` flag for the copies that are
recoverable only from multimappers.

## Usage

```
RUSTLE_VG_RESCUE_REPORT=1 rustle -L reads.bam --vg --genome-fasta ref.fa -o out.gtf
# → writes out.rescue.tsv
```

Opt-in (default-off → byte-identical when unset). Thresholds: `RUSTLE_VG_RESCUE_MIN_DECISIVE`
(default 3), `RUSTLE_VG_RESCUE_MIN_FRAC` (default 0.5).

## Columns (one row per VG-family transcript)

`chrom, tx_start, tx_end, strand, n_exons, coverage` — the isoform; `family_id, copy_id,
family_size, copy_locus` — the paralog copy; then the **decisive-evidence certificate**:

| column | meaning |
|---|---|
| `primary_reads` | reads whose PRIMARY alignment lands here — what a primary-only assembler sees |
| `decisive_own_reads` | `unique + strict` — reads that decisively belong to this copy (incl. owning multimappers) |
| `unique_reads` / `strict_reads` | place only here / multimappers that strictly own here by edit distance |
| `tied_reads` | multimappers tied with a sibling (non-identifiable mass) |
| `decisive_frac` | `own / (own + tied)` — 1.0 = unambiguous; low = phantom (e.g. DAZ3 ≈ 0.07) |
| `copy_independent_support` | the existing certified support fraction |
| `rescued` | recoverable only from multimappers (see classifier below) |

## The `rescued` classifier (TDD'd: `is_multimapper_rescued`)

A copy is rescued when it has **real, identifiable, multimapper-driven** evidence:

```
own_ev >= min_decisive          (a real copy, not noise)
AND decisive_frac >= min_frac   (identifiable — not a phantom whose reads tie to a sibling)
AND n_strict > n_primary         (strict-owning multimappers strictly EXCEED what a primary-only
                                  assembler sees — the rescue is real, not a copy StringTie sees)
```

The strict inequality matters: `own_ev = n_unique + n_strict`, and `n_unique` counts single-locus
*primary* reads already reflected in `n_primary`, so an `n_strict == n_primary` copy with unique
primaries is in fact primary-visible — not a multimapper-only rescue (caught in adversarial
review). Unit tests pin the corners: RABL2A-like (28 primary / 144 strict → rescued), invisible-real
(0 primary / 12 strict → rescued), phantom-tied (3 strict / 38 tied → NOT rescued), well-covered
(100 primary / 5 strict → NOT rescued), and the boundary (n_strict == n_primary → NOT rescued).

## Real-data demonstration (RABL2 family)

On the real RABL2 BAM (`--genome-fasta GGO.fasta`), the family forms across 3 copies and the
report flags **one rescued copy**:

| copy | primary_reads | decisive_own_reads | strict | tied | frac | rescued |
|---|---:|---:|---:|---:|---:|---|
| NC_073235.2 (15 isoforms) | 3 | 4 | 2 | 0 | 1.00 | false |
| **NC_086018.1** | 28 | **145** | 144 | 3 | 0.98 | **true** |

The rescued copy is decisively multimapper-owned (145 owners vs 28 primaries, frac 0.98 =
identifiable). Just as usefully, the report **exposes** that the 15 isoforms on NC_073235.2 rest
on only 4 decisive own-reads — a low-confidence signal a researcher should see. (Note: this is
rustle's *internal* `compute_copy_ownership` certificate over multimapping reads with a dNM
margin; the numbers are its own decisive-ownership metric, not the external strict-NM count from
the RABL2 worked example.)

## Implementation

- `vg.rs::CopyCertificate` + `is_multimapper_rescued` — pure, 5 unit tests.
- `pipeline.rs` — a read-only certificate block (`compute_copy_ownership` + `n_primary`) computed
  **before** the decisive gate mutates bundles; the report writer after `write_gtf`.
- Deterministic: rows sorted by `(chrom, start, family, copy)` (the GTF line-order's pre-existing
  rayon nondeterminism does not affect the report).

## Adversarial review (2-agent workflow)

Both lenses returned *ship-with-minor-fixes*; all fixed: **(major)** the classifier double-counted
unique primaries (`n_strict >= max(1,n_primary)` could mislabel a primary-visible copy at the
`n_strict == n_primary` boundary) → tightened to `n_strict > n_primary` + a regression test;
**(minor)** the row sort wasn't a total order (isoforms of one copy sharing a start tied, borrowing
the rayon iteration order) → full-row tiebreaker; **(minor/nit)** the rescued tally counted
transcript rows via a string-suffix probe → now counts distinct `(family, copy)` via `is_rescued()`.

## Verification

| check | result |
|---|---|
| unit tests | RED→GREEN, 6/6; full lib suite 268/0 |
| real RABL2 | 1 rescued copy flagged (144 strict > 28 primary, frac 0.98), log "rescued copies: 1" |
| determinism (2 runs) | `rescue.tsv` identical (total-order sort) |
| default (no flag) | no report file; GTF byte-identical |
