# Design: gated read-coherence layer (`--read-coherence`)

**Date:** 2026-06-15. **Status:** design for review.
**One-liner:** Productize read-coherent (per-molecule) path extraction as a gated, default-off
"better decision" layer that recovers real isoforms StringTie misses — **+580 exact (FSM) +
2,784 broad (=/c/j) annotated transcripts genome-wide** (validated) — by adding read-chains over
the flow floor and filtering the ~half-noise with an **annotation-free realness gate**
(canonical junctions + RT-switch + read-depth) plus **degradation-aware collapse**. Byte-identical
when off; `RUSTLE_PRECISE`-exempt.

## 1. Why (validated)
StringTie/flow assembly merges or drops real isoforms that individual long reads span. Read-chain
extraction (`RUSTLE_READCHAIN`) recovers them: rigorous cut (rc vs RefSeq minus StringTie vs
RefSeq, 25 chroms) = **+580 FSM** exact-annotated transcripts StringTie misses + **2,784** broad.
FSM ⇒ canonical & real by construction (not the PSV non-canonical-artifact failure mode). This is
~100–200× the PSV margin. **Cost:** the 580 sit inside ~10,144 raw read-chain extras (~half noise:
non-canonical, RT-switch, ISM-degradation), so a precision gate is mandatory to ship it.

## 2. Existing building blocks (reuse)
- `global_flow.rs::extract_transcripts_readchain` (727) — per-molecule chain extraction (source
  `"readchain"`). The generation engine. Today additive over `-G` flow (0 flow chains lost).
- `global_flow.rs::compute_5p_degrade_folds` (629) — folds 5′-truncated chains into their
  full-length parent (`RUSTLE_READCHAIN_DEGRADE_TOL`, default-on). **Incomplete:** 5′ only.
- `genome.rs` GenomeIndex — FASTA access (for canonical / RT-switch, both annotation-free).
- Coverage/`longcov` on transcripts — read-depth for the gate.
- **Absent (must add):** canonical-junction (GT-AG/CT-AC) check; RT-switch (direct-repeat) check.

## 3. Architecture (3 parts, all gated under `--read-coherence`)

### (a) Degradation-aware collapse (generation fix — kills the biggest noise bucket at source)
Extend `compute_5p_degrade_folds` → `compute_degrade_folds`: also fold **3′-truncated** and
**internal-fragment** chains into their full-length parent (a chain whose junction set is a
contiguous sub-path of a longer chain, sharing the appropriate terminus within `tol`), summing
abundance. Removes ~2,128 ISM-degradation extras *before* they are emitted (vs filtering after).

### (b) Annotation-free realness gate (filter — removes data artifacts)
Applied to read-chain-sourced transcripts only (never touches flow/guide transcripts), before
emission:
- **Canonical junctions:** every junction is GT-AG or CT-AC (FASTA donor/acceptor dinucleotides).
  Drop transcripts with any non-canonical junction (≈494 noise). New `genome` helper.
- **RT-switch:** drop a junction flanked by a direct repeat (≥`R` bp identical on the exon side of
  the donor and the intron side of the acceptor — the template-switch signature; ≈1,538 noise).
  New `genome` helper.
- **Read-depth:** require `longcov ≥ K` (default 2) — drops single-read intergenic/antisense
  (≈half of ~1,331).
Net: ~10,144 raw → ~4,000–4,700 real (keeps the 580 FSM + the canonical NIC/NNC), drops the noise.

### (c) Additive merge (strictly over the flow floor)
Gated read-chain transcripts that pass (a)+(b) and whose intron chain is absent from the flow
output are appended (additive; never displace flow/guide). Tag `source "read_coherence"`. Runs
through the existing identical-transcript dedup + TPM.

## 4. Flag / gating / invariants
- **Flag:** `--read-coherence` (CLI) → `RUSTLE_READ_COHERENCE=1`. Composes with `-G`.
- Internally enables `RUSTLE_READCHAIN` + the degrade-collapse extension + the gate.
- **Default OFF → byte-identical** to current output (the whole layer is skipped).
- **`RUSTLE_PRECISE`-exempt**; **strictly additive** (flow floor unchanged).
- Raw `--read-chain` stays as-is (ungated, for diagnostics).

## 5. Components / files
- `src/rustle/global_flow.rs`: extend `compute_5p_degrade_folds`→`compute_degrade_folds` (3′+internal);
  add the realness gate over read-chain transcripts at the readchain emission site.
- `src/rustle/genome.rs`: `is_canonical_junction(chrom,donor,acceptor)` + `is_rt_switch(...)` (FASTA).
- `src/rustle/main.rs`, `types.rs`, `pipeline.rs`: `--read-coherence` flag plumbing (mirror existing).
- **Reused unchanged:** flow assembly, dedup, TPM, emission.

## 6. Validation
- SQANTI3 (structural categories + canonical + RT) on the gated output via the copy-recovery harness.
- The rigorous cut: gated-rc-FSM(=) minus StringTie-FSM(=), genome-wide. **Success: keep ≈the 580
  FSM** (and a large share of the 2,784 broad) while the raw 10,144 → ~4–5k (noise dropped), with
  precision (FP rate) acceptable vs the ungated `--read-chain`.
- Per-stage: each gate component's keep/drop counts logged; default-off byte-identity re-verified.

## 7. Testing (TDD)
- Unit: `is_canonical_junction` (GT-AG/CT-AC true; others false); `is_rt_switch` (direct-repeat true).
- Unit: `compute_degrade_folds` folds a 3′-fragment and an internal-fragment into the parent (5′ still works).
- Unit: gate drops a non-canonical / RT-switch / low-depth read-chain tx; keeps a canonical supported one.
- Invariant: `--read-coherence` off → byte-identical to current; `RUSTLE_PRECISE` unaffected.
- Fixture: a locus where one read spans a full isoform the flow merges → `--read-coherence` emits it.

## 8. Risks
- **Gate over-kills the 580:** mitigated — FSM are canonical + supported, so the gate keeps them by
  construction; validate the 580 survive each iteration.
- **Precision cost:** read-coherence is additive; the gate is the precision control. If FP still too
  high, tighten depth/locus-context (the ~1,331 intergenic/fusion are the residual).
- **RT-switch heuristic imperfect:** it's a known-imperfect signal; tune `R` against SQANTI3's RTS calls.
- **Degrade-collapse correctness:** must not fold genuinely-distinct short isoforms — require exact
  sub-path + terminus-within-tol (not mere overlap); covered by unit tests.
