# Mis-chain Read Salvage — Design

**Date:** 2026-07-23
**Status:** design (approved: Approach A, keep-both)

## Goal

Recover the local alignment of reads that minimap2 mis-chained across a large spurious intron, instead of
discarding them, so that a real duplication locus is not silently lost. This must be a **general** fix for any
multi-copy gene family — validated on synthetic ground-truth and known families, NOT tuned to the Soto benchmark.

## Background — the root cause (verified)

In near-identical segmental-duplication families, minimap2 sometimes chains one read's two ends across a giant
gap, producing an alignment whose CIGAR carries a huge `N` (a spurious intron) linking two loci. Today such a
read is lost at one of two points:

1. **`pass1_skeletons_robust`** rejects any skeleton whose span exceeds `MAX_SPLICED` (300 kb). A read whose
   intron-chain contains a > 300 kb bridge produces an over-long skeleton → dropped before assembly.
2. **`retain_non_mischain`** drops the assembled bridging *transcript* (giant intron, 50 kb–300 kb, carried by
   `< MISCHAIN_MIN_JUNCTION_READS` reads at that exact junction).

Either way the read's **local** alignment at the member's locus is lost, because its intron-chain carries the
giant bridge and it groups into the bridging skeleton rather than a local one. When the mis-chained read is the
family's *other* paralog, its loss orphans the clean survivor, which is then discarded by the ≥2-loci
family-definition requirement (confirmed: AC134878.2 with 485 clean reads → 0 copies in isolation; recovered when
paired with its paralog TEKT4P2). This mechanism accounts for the 12 mis-chain Soto members + 3 collateral
singleton-drops = 15 members, and generalizes to any family with a hard paralog.

## Approach A — split the read before `pass1` (chosen)

Before skeletons are built, transform the read set: any read whose intron-chain contains a **spurious giant
intron** is split at that junction into local reads. The local segments then group into local skeletons and seed
normally. The bridging skeleton/transcript never forms, so both drop points become moot.

**Spurious giant intron** uses the *same* criterion as the existing gate (no new magic numbers):
`(acceptor - donor) > MISCHAIN_GIANT_INTRON_BP` (50 kb) **AND** the exact junction `(chrom, donor, acceptor)` is
carried by `< MISCHAIN_MIN_JUNCTION_READS` (3) reads in the original read set. A well-supported large intron
(≥ 3 reads — a real large-gene intron) is **never** split.

**Keep both segments:** both flanks were aligned locally by minimap2 with real homology (that is why it chained
them), so both are valid local reads. Downstream identifiability gates (≥2-loci, PSV, exon-sum refine) filter any
genuinely spurious locus. Splitting is symmetric — no "dominant flank" heuristic.

### The split function

```
/// Split reads at spurious giant introns (mis-chain salvage). A read whose intron-chain contains an intron
/// `(d,a)` with `a-d > giant_bp` AND `support[(chrom,d,a)] < min_reads` is cut at that intron into local
/// sub-reads; the giant bridge is removed. Reads with no such intron pass through UNCHANGED (clone). Support
/// is measured on the ORIGINAL read set (before any split). Deterministic; preserves input order, sub-reads
/// emitted in 5'→3' order in place of their parent.
pub fn split_mischained_reads(
    reads: &[PrimaryRead],
    support: &HashMap<(String, u64, u64), usize>,
    giant_bp: u64,
    min_reads: usize,
) -> Vec<PrimaryRead>
```

Algorithm, per read:
1. Exon boundaries are reconstructed from `(ref_start, introns, ref_end)`: exon `0 = [ref_start, d0]`, exon
   `k = [a_{k-1}, d_k]`, last exon `= [a_last, ref_end]`.
2. Mark each intron `(d_k, a_k)` as a **cut** iff `a_k - d_k > giant_bp` AND `support[(chrom,d_k,a_k)] < min_reads`.
3. If no cuts → emit the read unchanged (exact clone).
4. Otherwise partition the exon list at every cut. Each contiguous run of exons → one sub-read: `ref_start` =
   first exon start, `ref_end` = last exon end, `introns` = the non-cut introns internal to the run.
5. A single-exon run is valid (single-exon local read). Runs are emitted 5'→3'.

Multi-giant-intron reads split into ≥ 3 segments naturally (a cut at every sub-threshold giant intron).

### Wiring

Insert immediately before each `pass1_skeletons_robust` call in the detect paths (the gw-catalog paths at
`denovo_pipeline.rs:2024/2223/2364`, `detect_and_assign` at `1446`, and `detect_families` at `193`), guarded by
config:

```
let reads = if cfg.mischain_salvage {
    let sup = read_junction_support(&reads);          // support on ORIGINAL reads
    split_mischained_reads(&reads, &sup, giant_bp, min_reads)
} else { reads };
let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
```

The split feeds **seeding** (pass1 → skeletons → transcripts). The conflict/placement path re-fetches reads from
the BAM independently (`reads_in_region` at 2070 etc.), so copy-assignment is unaffected — the salvage changes
*which loci seed*, not read counts in assignment. `retain_non_mischain` remains as an idempotent safety net (after
the split there is nothing left for it to drop).

### Config

- `DenovoConfig.mischain_salvage: bool`, default **false**. `from_env()` reads
  `RUSTLE_MISCHAIN_SALVAGE=1` → true. **OFF ⇒ byte-identical to current output** (the split call is never made).
- Thresholds are the existing gate values, already env-overridable
  (`RUSTLE_MISCHAIN_GIANT_BP`, `RUSTLE_MISCHAIN_MIN_READS`) via `env_num`.

## Testing — the anti-overfit harness (runs BEFORE measuring the 15 targets)

1. **Unit — positive:** a synthetic `PrimaryRead` with exons `[100,200],[210,300]` bridged to `[80_000,80_100]`
   by a giant intron `(300, 80_000)` carried by 1 read → `split_mischained_reads` returns two reads:
   `[100..300, introns [(200,210)]]` and `[80_000..80_100, no introns]`.
2. **Unit — negative (mandatory):** the SAME giant intron but carried by ≥ 3 reads (a real large-gene intron)
   → the read passes through unchanged. Also: a read with only sub-`giant_bp` introns → unchanged. This proves
   the split fires only on spurious mis-chains.
3. **Unit — multi-cut:** two sub-threshold giant introns → three local sub-reads.
4. **Unit — pass-through byte-identity:** with `mischain_salvage=false`, a read set that WOULD split is returned
   verbatim (guards the OFF path).
5. **Integration — planted family:** simulate a 2-locus family where reads at locus B spuriously mis-chain to a
   distant decoy (giant intron). With salvage OFF, locus B is lost (family = singleton → 0 copies). With salvage
   ON, B seeds and the 2-locus family is recovered. (New fixture under `bench/` / `tests/`.)
6. **Regression — known families:** GSTM, MAGEA, DAZ, RBMY, TSPY, PCDHB — copy counts unchanged with salvage ON
   vs the committed baseline (these have no spurious giant introns, so ON must equal OFF).
7. **No new FP:** genome-wide (or a representative multi-chrom slice) family/copy count with salvage ON must not
   exceed OFF by more than the intended recoveries; spot-check that added copies are real loci, not bridges.
8. **Soto:** the members we already pass stay passing (no regression); then measure recovery among the 15
   mis-chain + collateral targets. This is the LAST check, not the objective.

## Success criteria

- OFF path byte-identical (test 4 + a byte-diff on a known catalog).
- Negative test passes (real large introns never split).
- Planted family recovered ON, lost OFF (test 5).
- Zero regression on the 6 known families and on already-passing Soto members.
- No genome-wide FP increase beyond intended recoveries.
- Net: recovers a meaningful share of the 15 targets (report the exact number honestly; partial is fine).

## Non-goals

- The inverted-duplication / same-strand fix (ID_260) — a separate spec.
- OR11H13P and other single-expressed-locus members — not bugs; handled by reclassification / reference-homology
  fallback, out of scope here.
- Changing the ≥2-loci requirement itself (we fix the upstream cause, not the gate).
