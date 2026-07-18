# Augment-and-Linearize Validation for Reference-Absent Copies — Design Spec

**Status:** approved (brainstorming), ready for implementation plan
**Date:** 2026-07-18
**Builds on:** the O4 absent-copy admission (`absent_copy.rs`, `Admission` enum); copy-graph v1+v2 (shipped, `34a672a`).
**Objective:** a threshold-free, falsifiable validation that a candidate reference-absent copy is real — "add it to the reference; do the previously-ambiguous reads become uniquely placed (linearize)? and does a length/composition-matched decoy fail to do so?"

## Goal

For each candidate reference-absent copy admitted by the existing gates, emit an auditable **linearize certificate**: augment the reference with the candidate copy, re-align the ambiguous (MAPQ-0) read pool, and measure the fraction that becomes uniquely placed (`MAPQ>0`) — compared against a dinucleotide-shuffled **decoy** of the same candidate. A real copy linearizes the pool; a decoy does not. Report the certificate by default; make it an opt-in admission gate.

## Motivation

Rustle already has a collapse-detection gate (`collapse_gate.rs`) but ships it **default-OFF** (`collapse_gate.rs:17-33`): a fixed genome-wide ambiguity background (`GENOME_WIDE_EPS_AMB = 0.001313`, 5785/4.4M MAPQ-0 primaries) cannot distinguish a genuinely reference-absent copy from generic unresolvable paralogy — *both* pile up MAPQ-0 reads, and a real paralog *already in the reference* also linearizes reads. The current admission gate (`absent_copy.rs:153-158`, remap identity < 0.98) only checks the candidate is *divergent enough*, never that it *resolves* the ambiguity. The **decoy control is the principled fix for exactly this confound**: a length- and composition-matched shuffle of the candidate is a null that a real copy beats but generic ambiguity does not. This is the operational, tractable form of "overlay the RNA on the genome graph; the overhang the genome can't cover is the missing copy; folding it in makes the reads linear" — reduced to sequence-to-reference alignment (doable) plus a permutation control (threshold-free, Canzar-style).

## Design

### §1 — The certificate

Per candidate copy, a `LinearizeCertificate`:
- `n_pool`: number of MAPQ-0 (ambiguous) reads at the family that we test.
- `linearized_frac_real`: fraction of the pool whose PRIMARY alignment lands on the **candidate contig** with `MAPQ>0` (the read found its true home on the added copy). This is a tighter, cleaner control than "MAPQ>0 anywhere": a decoy contig attracts ~no reads, so its fraction is a proper null.
- `linearized_frac_decoy` (mean over N decoys) + the per-decoy distribution.
- `delta = linearized_frac_real - mean(linearized_frac_decoy)`.
- `perm_p = (#decoys with linearized_frac_decoy >= linearized_frac_real + 1) / (N + 1)`.
- `verdict`: `LINEARIZES` (perm_p < α, default α = 0.05), `NOT` (perm_p ≥ α), or `UNDETERMINED` (pool too small / alignment failure).

"Unique" = `MAPQ > 0` — Rustle's native definition (minimap2 sets MAPQ=0 *exactly* when the primary is a coin-flip; `copy_assign.rs:683`, `denovo_assemble.rs:252`).

### §2 — Augment + re-align (real arm)

- **Augmented reference** = the family's current copy consensuses (`DenovoTranscript.seq` for each admitted copy) **+ the candidate `t.seq`** (the overlay transcript from gate 4, `absent_copy.rs:143-151`), written as a temp multi-FASTA (reuse the `project_families_batch` writer + `TempFile` RAII, `genome_projection.rs:14-15,129-140`). Contigs are intronless spliced consensuses, and the reads are exonic RNA reads, so alignment is exon-to-exon (no splice).
- **Pool** = the `region` reads (`denovo_pipeline.rs:1504-1513`, `AlignedRead.seq` present) whose `region_mapq == 0`.
- **Align** the pool read sequences to the augmented reference via minimap2 in PAF mode, reusing the `remap_identity_minimap2` shell pattern (`absent_copy.rs:223-299`): `minimap2 -c --secondary=no -t 1 <aug.fa> <reads.fa>` with a `map-hifi`/`asm20`-style preset (NOT `-x splice`; the consensus contigs are intronless). Respect `RUSTLE_MINIMAP2`.
- **Parse** each read's primary hit: PAF col 5 (target contig), **col 11 (MAPQ — no existing parser reads it; add it)**, primary via `tp:A:P`. Reuse the `parse_paf_hits` split loop (`genome_projection.rs:55-77`) extended with `f[11]`. `linearized_frac_real = (#pool reads whose primary hit is the candidate contig with MAPQ>0) / n_pool`.

### §3 — Decoy control

- Generate **N dinucleotide-preserving shuffles** of the candidate `t.seq` (same length, same di-nucleotide frequencies, scrambled order — a small deterministic generator, ~15 lines; seed from a hash of `t.seq` so certificates are reproducible). Default **N = 19** (permutation-p granularity 0.05); configurable.
- Add the **reverse-complement** of `t.seq` (`seq_utils::reverse_complement`, `seq_utils.rs:9`) as one extra cheap sanity decoy.
- For each decoy, build a **decoy-augmented reference** = family copy consensuses + the decoy contig (same contig count and lengths as the real arm — a matched null), align the SAME pool, and compute `linearized_frac_decoy_i` = fraction whose primary hit is the *decoy* contig with `MAPQ>0` (expected ≈ 0, since reads don't match a shuffle). (The decoy-augmented reference is the null; no separate baseline needed.)
- `perm_p` per §1. A real copy: the pool matches the real contig → `linearized_frac_real` high, decoys ≈ 0 → Δ large, `perm_p` small. A spurious candidate: real ≈ decoys → Δ ≈ 0, `perm_p` large.

### §4 — Report-first, opt-in gate

- **Default (report-only):** compute the certificate for every admitted candidate; admission is UNCHANGED (the existing remap-identity gate still decides). Emit `<out>.linearize.tsv`: `family_id, candidate_locus, n_pool, linearized_frac_real, mean_decoy, delta, perm_p, verdict`.
- **Opt-in `--linearize-gate`:** admission additionally requires `verdict == LINEARIZES`; a candidate that passes the remap-identity gate but is `NOT`/`UNDETERMINED` becomes a `DnaNeeds` record (needs DNA validation). This lets the certificate be validated on known-real (DSFAM26 MHC) and known-decoy cases before it decides anything.

### §5 — Hermetic architecture (testable without minimap2)

Mirror the existing `remap_identity` injection (`absent_copy.rs:115-159`). The certificate core:
```rust
pub fn linearize_certificate(
    candidate_seq: &[u8],
    family_copy_seqs: &[Vec<u8>],
    pool_reads: &[Vec<u8>],        // MAPQ-0 read sequences
    n_decoys: usize,
    realign: impl Fn(&[Vec<u8>] /*ref contigs*/, &[Vec<u8>] /*reads*/) -> Vec<Option<(usize /*best contig idx*/, u32 /*mapq*/)>>,
) -> LinearizeCertificate
```
The core builds the real/decoy contig sets, calls `realign`, and computes the fractions + `perm_p` — pure, no I/O. Unit tests pass a synthetic `realign` closure. The production wiring passes a real closure that shells minimap2 (temp FASTA + PAF parse). The candidate contig is a KNOWN index in the ref set, so "linearized on the new copy" is checkable.

### §6 — Wiring seam

At the stage-2 admission call site (`denovo_pipeline.rs:1524-1540`), on the `Admission::Copy(t, id)` arm: build the pool = `region` filtered by `region_mapq == 0` (sequences via `region[i].seq`), the family copy consensuses = the current `all_copies` `.seq`, and call `linearize_certificate` with the real minimap2 realign closure. Store the certificate (push to a `Vec<LinearizeCertificate>` returned alongside the families) and, under `--linearize-gate`, convert a non-LINEARIZES verdict into a `DnaNeeds`. The `<out>.linearize.tsv` writer lives in `copy_assign.rs` alongside `dna_needs.tsv`.

### §7 — Error handling / determinism

- Pool < `min_pool` (default 5 MAPQ-0 reads) → `UNDETERMINED` (no verdict; report NA). Never a false LINEARIZES on thin evidence.
- minimap2 spawn/exit failure → `Ok(None)` per read → `UNDETERMINED` (Rustle's graceful-degrade contract, `genome_projection.rs:48`).
- The decoy shuffle uses a **deterministic** seed (hash of `candidate_seq`) so the certificate is reproducible run-to-run (avoids the HashMap-ordering non-determinism flagged elsewhere).
- Runs only under `--absent-copies` (where candidates exist); default `copy_assign` output byte-identical otherwise.

### §8 — Testing

- **Hermetic unit tests** (synthetic `realign` closure): (a) real copy — pool reads map uniquely to the candidate contig, decoys get none → Δ large, `perm_p` small, `LINEARIZES`; (b) null candidate — real ≈ decoys → Δ ≈ 0, `perm_p` large, `NOT`; (c) pool < min → `UNDETERMINED`.
- **Dinucleotide-shuffle test:** the shuffle preserves length and di-nucleotide frequencies while changing order; deterministic for a given seed.
- **`perm_p` test:** the formula `(#≥ + 1)/(N+1)` at boundary cases.
- **Data-gated real render:** on DSFAM26 (real MHC collapsed copy), the certificate should be `LINEARIZES` (real ≫ decoys); a manually-shuffled candidate should be `NOT`. Foreground, small region, `winloci_scratch` (WSL2 crash rule).

## Out of scope

- Making the certificate the SOLE/replacement admission criterion (kept as report-first + opt-in gate; promotion to default gate is a later decision after validation).
- Genome-scale re-alignment against the whole reference (we augment the *local family* reference — the ambiguity is within-family).
- Reviving `collapse_gate`'s default-OFF binomial test (the decoy control supersedes its confound; we don't need the fixed-background version).

## Global constraints

- Default `copy_assign` (no `--absent-copies`) byte-identical; the certificate runs only in the absent-copy path.
- The certificate core is a pure function with an injected realign closure (minimap2-free tests).
- "Unique/linearized" = `MAPQ > 0` (Rustle's native definition); primary reads only.
- Deterministic decoy shuffle (reproducible certificates).
- Report-first by default; admission decisions unchanged unless `--linearize-gate`.
- No new k-mer computation; reuse the existing minimap2 shell + PAF-parse infra.
