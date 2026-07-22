# Multi-mapper (tied-secondary) seeding — Phase 1 (benchmark-validated)

**Date:** 2026-07-22
**Status:** approved design, pre-implementation
**Type:** pipeline feature (opt-in), + benchmark evaluation.

## Motivation

The Soto member-detection metric marks 86 members "missed" (recall 76.2%). Measured coverage shows
**all 86 have reads at their locus — none are truly silent.** 13 are covered *almost entirely by
AS-tied secondary reads* with 0–2 primaries (e.g. PPIAL4A: **0 primary / 780 secondary**). Their copy
plainly exists; every read simply maps equally well to a sibling (MAPQ 0), so it is *primary at the
sibling, secondary here*. The de-novo step seeds loci from **primary reads only**
(`pass1_skeletons_robust`), so a starved locus gets 0 credit and is called "missed."

K=0 (reads tie) fundamentally blocks **assignment** (which read is this copy's), not **existence**
(is there a copy here). This feature recovers the *existence* by seeding a locus from the agreeing
tied secondaries, marked **unassignable**, so detection recall reflects what the reads actually show
while the assignment metric is untouched.

## Goal (Phase 1)

Seed a candidate locus from agreeing AS-tied secondary reads when no primary anchors it, **on the
Soto benchmark first** (reference-guided is low-risk: the loci are known), and **measure**:
- recall gain — tied-seeded loci overlapping a previously-missed member (TP),
- over-seeding — tied-seeded loci overlapping no member (FP).

If FP is low and recall climbs, the gate is proven and we promote to genome-wide de-novo (Phase 2,
separate spec). Off by default; when off, behavior is **byte-identical** to today.

## Method

### New: `tied_seed_skeletons` (`denovo_assemble.rs`)

```
pub fn tied_seed_skeletons(
    tied_reads: &[PrimaryRead],       // AS-tied secondaries (already gated by tied_secondary_reads)
    primary_skeletons: &[Skeleton],   // pass-1 output, to dedup against
    min_reads: u32,
) -> Vec<Skeleton>
```

1. Group `tied_reads` by `(chrom, intron-chain)` — the **exact same keying** primary seeding uses
   (`introns: Vec<(u64,u64)>`). This is the shared-intron-chain agreement gate: a scattered pile of
   secondaries with no common junction structure forms no group of size ≥ `min_reads`.
2. For each chain-group with `>= min_reads` reads, compute its extent (min start / max end, honoring
   the same `min_terminal_support` trim as primary seeding).
3. **Distinct-locus dedup:** drop a group whose span overlaps any `primary_skeletons` entry on the same
   chrom (that locus is already seeded — we only add *starved* loci). Overlap = spans intersect.
4. Emit the survivors as `Skeleton { …, tied_seeded: true }`. Deterministic (sorted by
   `(chrom, introns)`), mirroring `pass1_skeletons_robust`.

### Struct / config changes

- `Skeleton` gains `pub tied_seeded: bool` (primary skeletons set `false`). All existing
  `Skeleton { … }` literals updated to `tied_seeded: false`.
- `DenovoConfig` gains `pub tied_seed: bool` (default `false`).
- `detect_families` gains a `tied_reads: &[PrimaryRead]` parameter:
  ```
  pub fn detect_families(reads: &[PrimaryRead], tied_reads: &[PrimaryRead],
                         genome: &GenomeIndex, cfg: &DenovoConfig) -> DenovoResult
  ```
  When `cfg.tied_seed && !tied_reads.is_empty()`, append `tied_seed_skeletons(tied_reads, &skeletons,
  cfg.pass1_min_reads)` to `skeletons` before `assemble_gate`. Otherwise the slice is untouched →
  byte-identical. Callers that don't use the feature pass `&[]`.

### Caller wiring (`gw_family_catalog` / `copy_assign` de-novo path)

The driver already fetches primary reads; when `--tied-seed` is set it also fetches AS-tied
secondaries via the existing `tied_secondary_reads` / `tied_secondary_reads_in_region`
(`as_ratio` reuses the existing `--as-ratio`, default 0.98) and passes them as `tied_reads`. A new
CLI flag `--tied-seed` sets `cfg.tied_seed`.

### Output / marking

`tied_seeded` flows from `Skeleton` onto the assembled transcript and its member/locus record.
Detection counts a tied-seeded locus as **detected**; its status string is **"detected (tied),
unassignable."** The identifiability/PSV gate downstream is **unchanged**: a tied-seeded locus has no
copy-specific PSV, so copy-assignment abstains on it — it can never falsely acquire read assignments.

## Phantom guard (three stacked gates)

1. **AS-tie** (`tied_secondary_reads`, ratio 0.98) — only genuine co-located ties, homology-shadow
   spillover already excluded.
2. **Shared intron-chain** — the transcript structure must agree (≥ `min_reads` on one chain).
3. **Distinct locus** — must not overlap a primary skeleton (no double-count).

## Benchmark evaluation (Phase-1 success gate)

Rebuild the Soto de-novo detection with `--tied-seed` on the scoped Soto BAM (`soto_reads.bam`), re-run
`soto_detection_eval.py`, and emit `bench/soto/tied_seed_eval.tsv`:
- recall before/after (276/362 → ?),
- per tied-seeded locus: overlaps a Soto member (TP) or not (FP),
- aggregate TP / FP → the over-seeding number.

**Success:** recall rises (the covered-but-tied members recover) with FP small (target: FP ≤ the recall
gain, i.e. precision of the new loci > 50%; report the actual number regardless).

## Testing (TDD)

Unit tests on `tied_seed_skeletons`:
- an agreeing chain (≥ `min_reads` sharing one intron-chain, distinct locus) → seeds one skeleton with
  `tied_seeded == true`;
- scattered secondaries (no chain reaches `min_reads`) → seeds nothing;
- a chain-group whose span overlaps a primary skeleton → deduped away (0 seeded);
- below `min_reads` → nothing.
Plus a `detect_families` off-path test: `tied_seed == false` (or empty `tied_reads`) → identical
skeleton set as before. Then the Soto eval for the recall / FP numbers.

## Out of scope

- **Phase 2** genome-wide de-novo promotion — separate spec, gated on Phase-1 FP being low.
- Any change to copy-assignment / PSV / significance logic — tied-seeded loci stay unassignable by the
  existing gate; we do not attempt to assign their reads.
- Reclassifying the other 73 primary-bearing misses (seeding/structural) — a different lever.
