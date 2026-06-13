# EM Restriction + PSV-FASTA Emission — Design

**Date:** 2026-06-12
**Status:** Design (awaiting review → implementation plan)

## Goal

Two related changes to Rustle's VG (paralog-family) mode, both motivated by the finding that
the fingerprint-EM is *largely redundant with NM-argmax for read assignment* and only adds value on
a narrow class of paralog families:

- **Part A — Restrict the EM** to run only where it can change the answer (genuine multimapper
  ambiguity + an anchor to resolve it), abstaining to 1/NH elsewhere. New VG default; legacy full-EM
  behind an opt-out. Contingent on a measurement showing net-neutral-or-better.
- **Part B — PSV FASTA** (opt-in, default-off): emit a spliced-sequence FASTA of every emitted
  transcript plus a companion `.psv.tsv` listing the copy-distinguishing variants — proving, against
  sequence-blind StringTie, that real PSVs were found.

These are independent and can be implemented/landed separately, but share the PSV machinery
(`enumerate_diagnostic_sites`) and the same opt-in-report idioms.

---

## Background / evidence (from investigation wf_7a251b82)

`run_fingerprint_em` (`src/rustle/vg.rs:4838`) reweights each multimapping read's copy assignment
(`reads[].weight`, `vg.rs:5378`), which flows into node coverage (`map_reads.rs:331`) → flow
decomposition (runs at `pipeline.rs:~11355`, *before* assembly). With the anchor prior on (default)
the E-step runs exactly once (`effective_max_iter=1`, `vg.rs:5139`).

Measured on chr19 (`GGO_19.bam`, `--vg --vg-snp --genome-fasta GGO.fasta`, vs `GGO_19.gtf`):

| config | tx | struct Δ | Sn / Pr |
|---|---|---|---|
| EM on (baseline) | 1995 | 0 | 94.7 / 87.3 |
| `--vg-em-iter 0` (anchor forces 1 iter — **not** a real disable) | 1995 | 0 | 94.7 / 87.3 |
| `--vg-em-uniform-prior` (**dead flag, never read**) | 1995 | 0 | 94.7 / 87.3 |
| EM truly off (`RUSTLE_VG_ANCHOR_PRIOR=0 --vg-em-iter 0`) | 2002 | 9 | 94.8 / 87.1 |

The EM fires on only **2 families** on chr19. Fully disabling it adds 7 transcripts at one
mixed-strand locus (`NC_073243.2:94.2–94.38 Mb`): **+2 TPs (class `=`) and +5 FPs** (4 `u`, 1 `x`).
Its non-redundant value: (a) the abundance prior soft-apportioning **truly-tied** multimappers
(prevents the documented DAZ/RBMY over-counting, `vg.rs:4933-4939`), and (b) recovering true unequal
copy split when **≥3 PSVs** exist. For plain read assignment it is redundant with the NM-argmax that
already picks the best placement (`vg.rs:5065`); most multimappers cover **0 PSVs**. Its real value
is concentrated on chrY paralog families (RBMY/TSPY/DAZ), not chr19.

Removing the EM entirely also drops the `capacity_confidence`/`copy_confidence`/`abundance_min`
annotations (`gtf.rs:191/202/213`), the tandem weight-floor (`pipeline.rs:~11767`), and mosaic
detection (`RUSTLE_VG_MOSAIC_ON`, hosted inside the EM). → **restrict, do not remove.**

---

## Part A — EM restriction

### A1. The gate

Slot a gate into the per-family loop of `run_fingerprint_em`, immediately after the fingerprint is
built (`vg.rs:~4909-4944`, where `fp.n_sites` is in scope) and before the `PreEntry` collection
(`vg.rs:~4958`). When the gate fails, push `EmResult::default()` and `continue` — reads keep their
incoming 1/NH weights, i.e. **byte-identical to the EM never having run for that family**.

Gate (run the EM body only if):

```
(fp.n_sites >= K) || any_anchor)        // condition A: there is something to resolve WITH
  && total_tied >= 1                     // condition B: there is genuine ambiguity to resolve
```

- `fp.n_sites` = `enumerate_diagnostic_sites` PSV count (already computed per family).
- `any_anchor` = at least one copy with `n_strict >= 1` **or** `n_unique >= 1` from
  `compute_copy_ownership` (`vg.rs:1857`) — a decisive NM-owner or a uniquely-mapped read gives the
  abundance prior a non-degenerate anchor. (Including `n_unique` is deliberate: the chr19
  precision-favoring family has 0 PSVs; whether it keeps the EM depends on whether it has any anchor
  — this is the **measurement-sensitive** knob, see A4.)
- `total_tied` = sum of `n_tied` across copies. If every multimapper is decisively owned, argmax
  already decides them → EM is pure redundancy → abstain.

`K` is exposed via the existing `RUSTLE_VG_EM_MIN_DECISIVE_SITES` (default 2; the gate uses the same
value). `compute_copy_ownership` is already invoked by `anchored_mass_per_copy` (`vg.rs:5120`); reuse
its `n_unique/n_strict/n_tied` tallies rather than recomputing.

### A2. `--vg-em-iter 0` cleanup

Today `--vg-em-iter 0` is silently overridden to 1 when the anchor prior is on (`vg.rs:5139`). Add an
early `if max_iter == 0 { results.push(EmResult::default()); continue; }` at the top of the per-family
loop so the documented iteration knob actually disables the EM independent of the prior. (Also fix or
remove the dead `--vg-em-uniform-prior` flag — parsed but never read; out of scope to wire its
intended behavior, but it must not look functional. Minimal: document it as deprecated/no-op, or
delete the CLI arg. Decision: delete the dead arg to avoid a misleading lever.)

### A3. Flags / default

- **New default:** the gate is ON. VG runs the EM only where it qualifies.
- **Opt-out:** `RUSTLE_VG_EM_LEGACY=1` → skip the gate, run the full EM on every eligible family
  (exact current behavior) for A/B comparison and rollback.
- The escape hatch is unaffected: EM is VG-only; `RUSTLE_PRECISE`/precise-mode is non-VG.

### A4. Measurement gate (REQUIRED before the default flip is finalized)

Because making the gate default changes VG output, the implementation MUST measure on chr19 and only
keep the default-on if net-neutral-or-better:

1. Run baseline (legacy), gated-default, on `GGO_19.bam`; compare transcript count, renumber-invariant
   coord-signature struct diff, and gffcompare Sn/Pr vs `GGO_19.gtf`.
2. If the gate causes the chr19 mixed-strand family to abstain and that regresses precision
   (the +5FP/−2TP swing), tune condition A (e.g. require `n_unique>=1` to keep that family in scope)
   or accept the swing only if Pr does not drop below baseline.
3. Spot-check a chrY paralog family (`bench/multi_copy_eval/ggo_Y.bam`, RBMY/TSPY/DAZ) to confirm the
   gate KEEPS the EM running there (those families have PSVs/ties → must pass the gate) and the
   over-counting guard is preserved.
4. **Fallback:** if net-neutral cannot be achieved as default, ship the restriction as opt-in
   (`RUSTLE_VG_EM_RESTRICT=1`) instead and leave the default as legacy. Record the outcome.

### A5. Tests

- Unit: gate predicate is a small pure helper `em_family_qualifies(n_sites, n_unique, n_strict, n_tied, K) -> bool`
  with cases: PSV-rich+tie → true; 0-site+owner+tie → true; 0-site+no-anchor+tie → false;
  any-evidence+no-tie → false; PSV-rich+no-tie → false.
- Unit: `max_iter==0` early-out yields `EmResult::default()` (no reweight).
- Integration/regression: `RUSTLE_VG_EM_LEGACY=1` reproduces pre-change VG output byte-for-byte on a
  small VG fixture.

---

## Part B — PSV FASTA (opt-in, default-off)

### B1. Flag

`RUSTLE_VG_PSV_FASTA` (env, default unset). When unset: no new files, GTF byte-identical. When set:
emit `<out>.transcripts.fa` and `<out>.psv.tsv` after `write_gtf`.

### B2. Emission

New block right after `write_gtf` (`pipeline.rs:19691-19692`), mirroring the `RUSTLE_VG_RESCUE_REPORT`
block (`pipeline.rs:19761-19799`). At that point `all_transcripts: Vec<Transcript>` is final (each has
`chrom`, `strand`, `exons: Vec<(u64,u64)>`) and `vg_snp_genome: GenomeIndex` is in scope.

For each emitted transcript:
- Build the spliced sequence: concatenate `genome.fetch_sequence(chrom, s, e)` (`genome.rs:113`) over
  sorted exons; if strand `-`, `reverse_complement` (`vg.rs:5952`).
- FASTA header: `>{transcript_id} family_id=.. copy_id=.. n_psv=N` (transcript_id from the same
  `assign_gene_tx_numbers` mapping the GTF uses, `gtf.rs:38`). `family_id`/`copy_id` present only for
  family copies; `n_psv=0` otherwise.

### B3. PSV recompute + join

PSVs are not alive at emit (the `ExonFingerprints` table is dropped inside the EM loop). Behind the
flag, recompute over the **original** `vg_families` (`pipeline.rs:~11496`): for each family,
`build_family_graph(fam, &bundles, genome, ...)` then `enumerate_diagnostic_sites(&fg, n_copies)`,
reading `per_copy_site_refs[copy_id]` → `(ref_pos, this_base, sibling_bases)`. Key the result by
`(family_id, copy_id)`, which equals transcripts' `(vg_family_id, vg_copy_id)` by the documented
`(family_id,copy_id)==fam_pos` invariant (`pipeline.rs:11064-11066`). Gate the recompute behind the
flag so default runs pay nothing.

### B4. `.psv.tsv` columns

`transcript_id  chrom  genomic_pos(1-based)  transcript_pos(1-based)  this_base  sibling_bases  family_id  copy_id`

`transcript_pos` is computed by walking the transcript's sorted exons accumulating length until
`genomic_pos` falls inside an exon (skip PSVs that fall in introns — they are not in the spliced
sequence). This makes the proof "this base IN this emitted sequence differs from the sibling copy."

### B5. Scope

Emit **all** transcripts' sequences (complete companion to the GTF). PSV-annotate only family copies;
non-family / single-copy transcripts get `n_psv=0` and no `.psv.tsv` rows. (Families with 0 diagnostic
sites — e.g. the segdup case where 229/230 reads had 0 PSVs, or mixed-strand 0-site fallbacks —
correctly report `n_psv=0`.)

### B6. Tests

- Unit: spliced-sequence builder (multi-exon, `+` and `-` strand RC) on a tiny synthetic genome.
- Unit: genomic→transcript coordinate mapping (PSV inside exon maps correctly; PSV in intron is
  dropped).
- Integration: `RUSTLE_VG_PSV_FASTA=1` on `GGO_19.bam` produces a FASTA with one record per GTF
  transcript and a `.psv.tsv` whose rows join 1:1 to family transcripts; default-off run is
  byte-identical (no new files).
- Showcase validation: run on a chrY paralog family and confirm `n_psv>0` rows exist (the
  StringTie-can't-see-this evidence).

---

## Invariants

- **Part B default-off → byte-identical** (no new files, GTF unchanged). Verified on a fixture.
- **Escape hatch** `RUSTLE_PRECISE=1` byte-matches 4705ab1 (both parts are VG-only; precise mode is
  non-VG, untouched).
- **Part A default flip is measurement-gated** (A4): net-neutral-or-better on chr19 + chrY paralog
  check, else fall back to opt-in.
- **`RUSTLE_VG_EM_LEGACY=1`** reproduces pre-change VG output exactly.

## Risks

- The gate could abstain on the chr19 mixed-strand family and shave precision (Pr 87.3→87.1). Mitigated
  by the `n_unique` anchor condition + the A4 measurement/fallback.
- PSV recompute rebuilds family graphs at emit — only under the flag, so no default cost; bounded by
  family count (small).
- `reverse_complement` / `fetch_sequence` coordinate conventions (0-based half-open exons vs genome
  fetch) must match the existing per-copy-sequence fetch in `build_family_graph` — covered by the RC
  unit test.
