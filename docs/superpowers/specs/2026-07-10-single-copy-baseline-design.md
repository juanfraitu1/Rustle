# Single-copy baseline & formal transcript — Design

**Date:** 2026-07-10. **Substrate:** gorilla (GGO) HiFi Iso-Seq. **Scope:** narrow — formalize the transcript
object, emit single-copy loci as the copy-number **baseline**, and internalize λ_global. **NOT** an isoform
assembler.

## Motivation

`depth_cn = E_fam / λ_global` is the project's genome-free copy-number claim (`readonly_copy_number.rs`), and
`λ_global` is *defined* as "median read count over single-copy transcripts." But today that scalar is computed by
an **external Python helper** (`bench/rna_copy_number_depth.py::global_single_copy_anchor`) and enters the shipped
Rust binary as a hand-supplied `--lambda-global` number. **A headline thesis quantity is computed outside the
tool.** This design closes that provenance gap: the single-copy transcript set and λ_global are computed inside
the binary, once, genome-wide, and cached.

It also formalizes the transcript object and the unifying claim that a **locus is a copy family of size ≥ 1** —
single-copy is the degenerate χ(H)=1 case of the same operator that resolves multi-copy families.

## Non-goals (the scope boundary that keeps this thesis-aligned)

The thesis is **not an assembler** (StringTie/network-flow retired at F1 ~93.5). This design must not drift back:

- **No genome-wide isoform GTF.** The single-copy output is a **TABLE** (no exon records, not IGV-loadable), so
  it does not invite `gffcompare` against StringTie/IsoQuant/FLAIR. The isoform-level `--gtf` stays an opt-in
  side-output, unchanged.
- **No isoform-level quantification.** Distributing an assigned read to a specific isoform within a copy is the
  StringTie within-locus deconvolution problem — explicitly out of scope.
- **No reconciliation of the `--gtf` k=1 vs assignment k=2 divergence** (a known separate inconsistency; noted,
  not fixed here).
- The single-copy table is framed only as **copy-number calibration**, never transcript-reconstruction accuracy.

## The formal transcript definition (documentation, no new type)

A **transcript** is the existing `DenovoTranscript { chrom, start, end, n_reads, strand, introns, seq }`
(`family_detect.rs:63`). Its definition is the `assemble_gate` predicate (`denovo_assemble.rs:759`), stated
precisely and written into the spec + a module doc-comment as the canonical definition:

> A transcript is an exact-intron-chain cluster of primary reads whose **locus** — the connected component of the
> junction-incidence graph (`locus_support`) — carries ≥ `GATE_MIN_READS` (3) reads, whose junctions are all
> canonical and consistent-strand, and whose spliced length lies in `[MIN_SPLICED, MAX_SPLICED]` = `[100,
> 300000]` bp.

No new `Transcript` type or trait is introduced — the object exists; this design names and documents it.

## The unification (a claim, not code)

A locus is a copy family of size ≥ 1:
- **single-copy** = χ(H) = 1: zero PSVs, every read trivially assigned to the one copy. `chi_h` on a one-copy
  allele set already returns 1 (`readonly_copy_number.rs`).
- **multi-copy** = χ(H) ≥ 2, resolved by the significance gate.

The per-region `copy_assign` output is **unchanged** — it still emits only multi-copy families. The unification
is a conceptual/defense claim; single-copy loci are surfaced only by the genome-wide baseline pass below, so the
family tables are never flooded (measured genome-wide ratio: ~102,456 transcripts : ~425 family copies ≈ 241:1).

## Component 1 — genome-wide single-copy baseline pass (a FLAG on `gw_family_catalog`)

`gw_family_catalog --single-copy-baseline` (new flag on the existing binary — it already owns the per-chrom,
free-reads-between-chroms genome-wide traversal and the OOM discipline).

**Behaviour:**
1. Traverse the genome per-chromosome (reuse the `detect_conflict_catalog_genome_wide` traversal, which already
   loads primaries per contig and `drop(reads)` before the next), using the **assignment-path assembly**
   (`pass1_skeletons_robust(min_terminal_support)` + `assemble_gate` with `pool_locus_support` + span-aware
   collapse) — so single-copy transcripts are **consistent with the copies**, not the divergent `--gtf` path.
2. **Single-copy = a rep that no emitted family claims.** The traversal already produces all reps and the family
   set; a single-copy locus is a rep whose tid appears in no `ColocatedFamily`.
3. Emit `<out>.single_copy.tsv`: `chrom  start  end  strand  n_reads  n_exons  chi_h(=1)  n_psv(=0)`, one row per
   single-copy locus, sorted by (chrom, start).
4. Compute `λ_global = median(n_reads)` over the single-copy loci and write `<out>.lambda_global.tsv`: a header
   line + one data line `lambda_global  n_single_copy_loci` (the scalar plus its sample size, for provenance).

**Interface change (the one real plumbing item), memory-safe by construction:** at ~100k transcripts genome-wide,
holding full `DenovoTranscript` (with `seq`) would blow memory and break the per-chrom discipline. So the
single-copy set is **filtered per chromosome, to a lightweight record**, and the full reps are freed with their
chrom. Introduce:

```
pub struct SingleCopyLocus { pub chrom: String, pub start: u64, pub end: u64,
                             pub strand: char, pub n_reads: u32, pub n_exons: usize }
pub fn single_copy_loci(reps: &[DenovoTranscript], families: &[ColocatedFamily]) -> Vec<SingleCopyLocus>
// a rep whose tid appears in no family; carries no seq
pub fn lambda_global(loci: &[SingleCopyLocus]) -> Option<f64>   // median n_reads; None on empty
```

The genome-wide traversal calls `single_copy_loci(chrom_reps, chrom_families)` **per chromosome**, appends the
lightweight records to a genome-wide accumulator, then drops that chrom's full reps (existing `drop(reads)`
discipline extends to reps). `lambda_global` runs once over the accumulator. `detect_families`'s `DenovoResult`
(`denovo_pipeline.rs:111`) is the per-region counterpart and is not changed.

**λ definition detail:** the median is over single-copy loci that passed the gate (`n_reads ≥ GATE_MIN_READS`).
A locus counts as single-copy if it forms no multi-copy family in the data — which includes truly single-copy
genes and any whose paralog is unexpressed/collapsed. For a normalization floor this is the correct population
(the median expression of loci that behave as single-copy), and the caveat is documented.

## Component 2 — λ internalization in `copy_assign`

Add `copy_assign --lambda-file <path>` that reads the scalar from a `lambda_global.tsv`. `depth_cn` is unchanged
— it still consumes one `f64`. Precedence: explicit `--lambda-global <f64>` (manual override) > `--lambda-file` >
absent (`depth_cn = NaN`, emits `NA`, exactly as today). The external Python helper is thereby retired from the
copy-number path (kept only as a cross-check, not a dependency).

## Files

- **Modify** `src/bin/gw_family_catalog.rs`: `--single-copy-baseline` flag; emit `.single_copy.tsv` +
  `.lambda_global.tsv`.
- **Modify** `src/rustle/vg_family/denovo_pipeline.rs`: a genome-wide function that returns single-copy reps
  alongside families (the interface change); document the transcript definition on `assemble_gate` /
  `DenovoTranscript`.
- **Create** `src/rustle/vg_family/single_copy.rs`: `SingleCopyLocus`, `single_copy_loci(reps, families) ->
  Vec<SingleCopyLocus>` (lightweight, no seq), and `lambda_global(loci) -> Option<f64>` (median n_reads).
- **Modify** `src/bin/copy_assign.rs`: `--lambda-file`; wire into the existing `depth_cn` path.

## Testing

- `chi_h` on a one-copy allele set returns 1 (pin the χ(H)=1 boundary case).
- `single_copy_loci`: reps `{a,b,c}`, a family `{a,b}` ⇒ single-copy = `{c}`; the returned record carries no seq.
- `lambda_global`: median of `[10, 20, 30, 40]` = 25; median of `[10,20,30]` = 20; `[]` ⇒ `None`.
- `--lambda-file` parsing: a `lambda_global.tsv` with `lambda_global=25` is read and `depth_cn(50, 25) = 2.0`.
- Precedence: `--lambda-global 30 --lambda-file <25-file>` ⇒ 30 wins.
- A small genome fixture (one single-copy locus + one 2-copy family) run through the baseline flag emits a
  1-row `single_copy.tsv` and the correct λ; the 2-copy family does NOT appear in `single_copy.tsv`.

## Reproduce

```
gw_family_catalog --bam GGO_mm.bam --fasta GGO.fasta --single-copy-baseline --out gw
# -> gw.single_copy.tsv  +  gw.lambda_global.tsv
copy_assign --bam ... --fasta ... --region <r> --lambda-file gw.lambda_global.tsv --out fam
# depth_cn now uses the in-binary genome-wide λ, no --lambda-global scalar, no python helper
```

Related: `bench/FAMILY_ARTIFACT_AUDIT.md` (the 241:1 scale), `readonly_copy_number.rs` (depth_cn/λ), the
retired StringTie parity era (why this stays a table, not a catalog).
