# Tied-secondary seeding — Phase-1 benchmark result (2026-07-22)

**Verdict: the feature is correctly implemented and safe, but it does NOT recover Soto members on this
benchmark.** An honest negative result, with the mechanism fully diagnosed below.

## What was built (Tasks 1–4, all reviewed clean, unit-tested)

`tied_seed_skeletons` (`denovo_assemble.rs`) groups AS-tied secondary reads by `(chrom, intron-chain)`,
keeps chains with ≥3 reads whose span doesn't overlap a primary skeleton, and returns them as
`Skeleton { tied_seeded: true }`. Wired into `detect_and_assign` behind `cfg.tied_seed` (CLI
`--tied-seed`); off = byte-identical. 4 unit tests pass; the secondary→`PrimaryRead` builder does
populate introns (so the intron-chain gate is not a silent no-op).

## Benchmark (direct `copy_assign` on Soto family regions)

Ran `copy_assign --recover-copies [--tied-seed]` on co-located Soto family regions and diffed the
detected copies (`families.tsv` `n_copies`). Result: **no change** — with the default 0.98 AS-tie gate
AND with a loosened 0.85 gate.

- TRIM64 region (chr11:89.68–90.00 Mb, 4 co-located copies, 2 missed): OFF = 13 copies; `--tied-seed`
  (0.98) = 13; `--tied-seed --as-ratio 0.85` = 13. No copy added.

## Why (two structural reasons, both correct behavior)

1. **The genuine-K=0 tied-starved members it targets are DISPERSED.** The covered-but-tied members with
   0 primary + heavy *genuine* ties (PPIAL4A: 0 primary / 780 secondary; DDX11L16; WASH6P) sit in
   families whose copies are 20–30 Mb apart. `copy_assign` detects **co-located** families within a
   ~5 Mb window, so it never assembles these into a family — a lone tied-seeded locus is not emitted.
2. **The co-located missed members are distinguishable-but-undercovered, not tied.** TRIM64B/TRIM64,
   DEFB108C/D, AC246785.1 have primary reads and were "recovered by coverage" in the top-up experiment
   (they have distinguishing sequence). Their secondaries mostly fail the 0.98 AS-tie gate. And where
   secondaries are admitted (0.85), they **don't agree on a spliced chain**: at TRIM64B the strongest
   agreement is the *unspliced* chain (13 reads — correctly skipped by the shared-intron-chain gate),
   while the spliced reads scatter across 70 distinct chains (max 5 sharing one, mostly 3). So the
   agreement gate correctly declines to seed a robust locus.

**Over-seeding: 0.** The feature added no spurious loci either — the phantom guards (AS-tie + shared
intron-chain + distinct-locus) hold.

## Interpretation

tied-seed is the right tool for **co-located, genuine-K=0** starved copies. Soto's missed members are
almost entirely **either** dispersed-K=0 (need a reference-guided / cross-locus detection path, not the
co-located `detect_and_assign` path) **or** distinguishable-but-undercovered (the coverage top-up
experiment recovered 12 of those — `soto_topup_recovery.tsv`). Neither niche is what tied-seed serves,
so on Soto it is correctly a safe no-op.

The feature ships opt-in (byte-identical off) and is ready for the co-located K=0 case and for a Phase-2
`gw_family_catalog` integration; it does not, on its own, move the Soto member-detection recall.

Reproduce: `target/release/copy_assign --bam A119b.t2t.bam --fasta chm13v2.0.fa --region
chr11:89676664-90003077 --recover-copies --tied-seed --skip-poa-diagnostic --out <scratch>` vs the same
without `--tied-seed`; diff `<out>.families.tsv`.
