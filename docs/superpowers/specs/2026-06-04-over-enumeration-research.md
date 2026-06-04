# VG over-enumeration: root cause + a novelty-safe suppression signal (research start)

**Date:** 2026-06-04. **Status:** research grounded; candidate signal validated on hidden-copy ground
truth; NOT implemented (pre-ship validations listed). 2-agent workflow `over-enumeration-research`.

## The problem
On the RBMY tandem mega-bundle, VG+tandem emits 23 transcripts (Tx Sn/Pr 55.0/47.8) vs baseline's 16
(65.0/81.2) — ~7–10 spurious extras kill precision (c4 → 9 transcripts vs ~7 real; c0 → 4 for 1 real).
A coverage filter cannot remove them (entangled with real low-cov copies). Goal: suppress the spurious
extras **without** killing genuine novelty (copies absent from the annotation; hard-to-model
reshufflings = inversions / gene-conversion).

## Pollution audit (the "don't let previous attempts confound it" requirement)
- **completion** (`RUSTLE_VG_COMPLETION_OFF=1`) and **borrow** (`RUSTLE_VG_NO_BORROW=1`): **INERT** — no
  change to the over-enum. Clean.
- **EM anchor-prior** (`RUSTLE_VG_ANCHOR_PRIOR=0 RUSTLE_VG_JOINT_STRAND_EM=0`): a genuine **confound** —
  turning it off IMPROVES to 65.0/54.2 (the sink-collapse that apportions 96/110 reads onto c4 bloats
  it by ~1 transcript). A control variable, not an assumption.
- The ~90 other experimental flags (hidden-copy, mosaic, segdup, inversion, FAMILY_BOOST, error-aware)
  are all default-OFF (verified) — not silently in the mix.

## Root cause: SECONDARY-READ CONTAMINATION (not novel structure)
tspy.bam has 87 primary vs **512 secondary** alignments. Cross-mapped secondaries (e.g. c0-origin
reads, mapq 60, landing on c4) draw splice junctions that **no primary read supports**, and
`path_extract` enumerates them into extra transcripts.
- gffcompare: **0/72 novel introns** — every extra is a *real* intron in a *contaminated combination*,
  not invented structure.
- The 8 phantom extras each carry a **secondary-only junction** (0 primary reads); **0/47** real
  (`=`-matched) junctions are secondary-only. Clean separation.
- Proof: `samtools view -F 0x100` (primary-only) → VG+tandem reproduces baseline **exactly** (16 tx,
  65.0/81.2). Every spurious extra vanishes with secondaries removed. EM sink-collapse is a minor
  (+1 tx) secondary contributor; `path_extract` is the downstream amplifier, not the root.

## The novelty-safe signal: own-chain compatible PRIMARY-read support (`compatPrim`)
`compatPrim(tx)` = number of **primary** reads whose intron chain is compatible (subset, ±6 bp) with
the transcript's chain — reads that actually carry *this* isoform's junctions.
- **Separation AUC 0.992** (vs coverage 0.777, vs raw locus primary-support 0.697 — the latter has no
  *within-copy* discrimination, which is where the over-enum lives). Spurious extras have
  `compatPrim = 0`; genuine isoforms have 3–16. Threshold ≥3 → all 11 genuine kept + 1 borderline-real
  → **Pr 91.7% (above baseline's 81.2%)**.
- **Novelty preserved (decisive test):** a synthetic **hidden copy absent from the reference**
  (`gen_synthetic.py --hidden-copies 2 --distinct-isoforms`) is surfaced with `compatPrim = 29` (its
  mismapped own-reads back its chain) and is **preserved at every threshold** → false-suppression of
  genuine novelty = **0**. Mechanism: a real novel copy emits real primary reads that form its chain;
  contamination artifacts do not. It uses **no annotation**, so it works for copies-absent-from-annotation.

## Honest scope limits (measured, not hand-waved)
1. **Gene-conversion / mosaic recombinants are BLIND to this signal** — they differ by *exonic SNPs*,
   not junctions, so they carry no novel junction and get pooled into the parent copy. These must be
   handled by the **PSV/SNP-phasing channel** (the diagnostic-site fingerprint), NOT this filter, or
   genuine recombinants are silently merged. → **Two-channel framework:** `compatPrim` for
   junction-level novelty (hidden copies, inversions, exon-skips); PSV channel for SNP-level novelty
   (gene conversion). They are orthogonal.
2. **Low-expression novel copies are lost UPSTREAM** (VG drops them at assembly before any filter).
   `compatPrim` is safe (won't mis-suppress) but cannot rescue what VG already dropped.
3. On real RBMY, `compatPrim≥3` drops all c0/c6 transcripts — but those VG assemblies were all *wrong*
   (FP anyway), so no annotated isoform is lost; locus recall on c0/c6 falls. → implement as **SOFT
   demotion** (flag/downweight, not a hard cut), with a guard exempting the **dominant primary chain of
   an otherwise transcript-less locus** so a real-but-mis-assembled copy still yields one transcript.

## Proposed direction (NOT yet implemented)
1. **Junction-admission gate in tandem `path_extract`:** require ≥1 primary read to admit a junction
   into the enumeration graph (or demote 100%-secondary junctions). Removes exactly the 8 phantom
   extras; strictly better than the blunt `-F 0x100` strip (which also discards secondaries that
   legitimately *confirm* primary-supported junctions).
2. **`compatPrim` soft-demotion** of emitted transcripts + the dominant-chain exemption guard.
3. Keep the EM anchor-prior **off for tandem** (it bloats the sink).
**Pre-ship validations:** ✅ **inversion DONE** — an `--invert-copies` inverted copy emits 3
transcripts (− strand) with `compatPrim` **40–47** (≫ threshold 3): its reads align as primaries at the
inverted locus and back their own chain, so the gate **preserves inversions** (same mechanism as the
hidden copy). Remaining: ≥2 more synthetic seeds; confirm gene-conversion stays on the PSV channel
(don't let `compatPrim` touch it). Do NOT promote c6 as recovery — its VG transcript is class
`n` (wrong structure, 2 primary / 62 secondary reads = the DAZ3 false-positive pattern).
