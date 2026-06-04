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

## Prototype BUILT + validated (2026-06-04, opt-in `RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS`)
Architectural finding from the first attempt: filtering `junction_stats` is **bypassed** — VG
assembly (`path_extract`) is **read-driven** (it re-derives junctions from the reads). So the gate
must act on the READS. Implemented (pipeline.rs, after the EM has used the secondaries for
apportionment + `copy_confidence`, scoped to `TandemCopy` sub-bundles): keep all primary reads + only
the secondary reads whose **entire chain is primary-supported** (they confirm real junctions → real
coverage); drop secondary reads carrying any **phantom (zero-primary) junction**.

| RBMY mega-bundle | tx | Tx Sn/Pr | Locus Sn/Pr |
|---|---|---|---|
| baseline (non-vg) | 16 | 65.0 / 81.2 | 83.3 / 100 |
| gate OFF | 23 | 55.0 / 47.8 | 66.7 / 66.7 |
| strip-ALL-secondaries (too blunt) | 8 | 40.0 / 100 | 50.0 / 100 |
| **gate ON (targeted)** | **15** | **55.0 / 73.3** | 66.7 / 66.7 |

Removes exactly the 8 phantom-junction transcripts → **precision 47.8 → 73.3 with ZERO recall loss**
(Sn stays 55.0). The blunt "strip all secondaries" over-corrects (Sn 55→40) because it discards the
coverage real low-cov copies need — the targeted version (drop only phantom-junction secondaries)
keeps it. **Novelty preserved** (the decisive test): a hidden copy absent from the reference (3 tx →
3 tx) and an inverted copy (3 tx → 3 tx) both survive the gate unchanged. Regression-safe: opt-in +
scoped to TandemCopy → gate-OFF unchanged (RBMY 23 tx), DAZ3 4.04, default headline 95.6/90.5, suite
222/0, obj5 1.0/abstain, tier1 100.

Remaining gap to baseline (Sn 55 vs 65) is NOT the phantom over-enum — it's c0 mis-assembly + the c6
class-`n` false positive (separate issues; do NOT promote c6). Next: ≥2 more seeds + a gene-conversion
case (must stay on the PSV channel, untouched by this junction-level gate) before any default change.

## c0 mis-assembly investigated (2026-06-04) — diagnosed, NOT cleanly fixable by the gate
After the gate removes the phantom over-enum, the dominant residual is c0 (LOC129530243).
Diagnosis:
- **c0 IS genuinely expressed and matches RefSeq**: 4–6 full-length primary reads carry its 9-intron
  chain, identical to XM_055377108 within ±1 bp. So it is a real mis-assembly, not annotation noise.
- **The failure is a RETAINED INTRON over 4 consecutive MICRO-EXONS** (89 / 111 / 111 / 115 bp): VG
  merges 19602756–19606688 into one exon instead of splicing the micro-exons.
- **It is tandem-specific, NOT a general micro-exon limit**: baseline (primary-only) SPLICES the
  micro-exons correctly. The tandem sub-bundle retains them because **secondary reads from sibling
  copies span the micro-exon region UNSPLICED** (their copy lacks the micro-exons) → retained-intron
  coverage that beats c0's 4–6 spliced primary reads → path_extract retains.
- The EM anchor-prior makes it worse (sink-collapse starves c0's primaries); **gate + EM-off reaches
  baseline transcript-Sn parity 65.0** overall, but c0 specifically still isn't recovered.

Attempted fix (drop secondaries that RETAIN a primary-supported intron): it **deletes c0 entirely**
(over-removal: c0 → 0 tx) rather than recovering it, and achieved "65.0/81.2 = baseline parity" only
by dropping the hard copy + over-enumerating c2 — a misleading aggregate. **Reverted** (dropping a
real expressed copy conflicts with the find-copies goal). **No config — including baseline — recovers
c0's exact chain.** The real lever is to make the tandem sub-bundle assembly emit c0's full-length-read
chain directly (as baseline's primary-only bundle does) — a deeper tandem-assembly / path_extract
question (why the sub-bundle assembly diverges from baseline's), out of scope for the junction gate.

## Proposed direction (original — now prototyped above)
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
