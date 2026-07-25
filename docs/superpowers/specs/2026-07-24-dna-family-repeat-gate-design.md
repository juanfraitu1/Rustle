# DNA-family-fallback Repeat/Complexity Gate — Design

**Date:** 2026-07-24
**Status:** design (approved: repeat gate on the orphan locus, soft-mask threshold)

## Goal

Make `dna_family_fallback` (the extended single-locus definition) honest by excluding re-admissions whose
re-admitted locus is repeat/transposon-derived — the failure mode an adversarial audit exposed (11/12
off-benchmark and both on-member controls were Alu SINEs / low-complexity, not real gene-family SDs). This
mirrors Soto's own definition, which explicitly RepeatMasks before calling an SD.

## Background — the verified problem

Genome-wide, `dna_family_fallback` re-admitted 715 loci (89.2% → 94.8% "recall"). An adversarial workflow +
deterministic soft-mask check showed re-admissions are ~57-59% RepeatMasker-soft-masked (only 11% complex-unique),
and of 10 members "recovered", only 1 (TRIM64, 29% soft-mask) survives a repeat filter. The current admission
gate (`≥1 projection at min_identity=0.90`, no repeat filter) projects orphan loci onto dispersed Alus (~99%
identical genome-wide), which masquerade as high-identity paralogs. The +20 recall is a repeat artifact.

## The gate

In `readmit_dna_family_batch` (`collapse_enumerate.rs`), before admitting a candidate, compute the soft-mask
fraction of its re-admitted locus and reject if repeat-dominated.

- **Read soft-mask verbatim:** use the existing `repeat_catalog::IndexedFasta` (reads bytes lowercase-preserved,
  purpose-built for soft-mask; `GenomeIndex` upper-cases and cannot be used). Open once per batch.
- **Signal:** `softmask_frac(seq) = lowercase_bases / total_bases` over the candidate's `(chrom, lo, hi)`.
- **Reject** the candidate when `softmask_frac >= max_softmask`.
- **Threshold:** `max_softmask` from `env_num("RUSTLE_DNA_FAMILY_MAX_SOFTMASK", 0.50)` — reject
  majority-repeat loci. 0.50 keeps the complex-unique real recoveries (TRIM64 @0.29) and drops the ~57-59% repeat
  FPs. Env-tunable so the recall/precision trade-off can be swept without recompiling.
- **Scope:** gate the ORPHAN locus (the object being re-admitted); the audit showed the orphans themselves are
  repeat-dominated, so gating them suffices. Projections are not separately gated (keeps it simple; can add later).
- **Inert unless active:** only runs inside `readmit_dna_family_batch`, which only runs when
  `cfg.dna_family_fallback` is set. Byte-identical for every default (fallback-off) run.

### Signature change

`readmit_dna_family_batch` gains the threshold parameter (threaded from `DenovoConfig`):
```
pub fn readmit_dna_family_batch(
    candidates: &[(String, String, u64, u64, Vec<u8>)],
    bam_path: &str, fasta_path: &str, minimap2: &str, threads: usize,
    min_identity: f64, max_softmask: f64,   // NEW
) -> Vec<ExpressedCollapsedFamily>
```
Add `pub fn softmask_frac(seq: &[u8]) -> f64` (pure, unit-tested).
Add `DenovoConfig.dna_family_max_softmask: f64` (default 0.50, `from_env` reads `RUSTLE_DNA_FAMILY_MAX_SOFTMASK`),
passed at the call site (`denovo_pipeline.rs:2470`).

## Testing

1. **Unit — `softmask_frac`:** all-lowercase → 1.0; all-uppercase → 0.0; half/half → 0.5; empty → 0.0; N/n handled.
2. **Unit — gate decision:** a candidate whose region is `>= max_softmask` lowercase is dropped; a complex
   (mostly uppercase) candidate with a valid projection is kept. Use a tiny in-repo test FASTA + `.fai` (or a
   pure-function gate helper `passes_softmask_gate(frac, thresh)` tested directly to avoid FASTA I/O in the unit).
3. **Byte-identity:** with `dna_family_fallback` off, output unchanged (the gate code is never reached).
4. **Re-measure (resumable per-chrom):** `dna_fallback_perchrom.sh` with the gate ON:
   - off-benchmark re-admission count must drop sharply (target: the ~508 → a small fraction),
   - the complex-unique real recoveries (TRIM64 etc.) survive,
   - report gated recall (expected ~89-90% honest, not 94.8%).

## Success criteria

- `softmask_frac` unit tests pass; gate drops repeat loci, keeps complex ones.
- Default (fallback-off) byte-identical.
- Gated per-chrom re-measure: off-benchmark volume collapses; report the honest recall + the real (complex-unique)
  members recovered. A small honest gain is the expected, acceptable outcome.

## Non-goals

- Gating projections (only the orphan locus, for now).
- Recovering the compound-failure members (ANKRD36B/CDH12/NCF1) — established as an RNA identifiability limit.
- Wiring `dna_family_fallback` into the default `--cross-chrom` Soto path (separate decision, after the gate proves honest).
