# RBMY1 — copy-level assignment validation on a testis Y gene family

**Question:** does the copy-assignment method actually work on a multi-copy
testis gene family — i.e. does it recover multiple expressed copies, or collapse
them? RBMY1 (RNA-binding motif protein, Y chromosome, family 1) is an ideal,
biologically-known test: a testis-specific Y ampliconic family.

## The locus

A 6-copy tandem array on chrY, `NC_073248.2:19,602,754–19,730,926` (~14 kb per
copy, − strand, ~23 kb spacing), annotated LOC129530238–243 with descriptions
"RNA-binding motif protein, Y chromosome, family 1 member B-like / F-J-like".

| label | gene | position |
|---|---|---|
| c1 | LOC129530243 | 19,602,754–19,616,644 |
| c2 | LOC129530239 | 19,625,715–19,639,621 |
| c3 | LOC129530240 | 19,648,670–19,662,577 |
| c4 | LOC129530238 | 19,671,638–19,685,525 |
| c5 | LOC129530241 | 19,694,606–19,708,531 |
| c6 | LOC129530242 | 19,717,578–19,730,926 |

**Pairwise identity:** the core c2–c5 are ~99.8% identical; c1 ~97.3%; c6 the
outlier at ~93–96%. A genuine identity gradient, as expected for a young amplicon.

## Reads, multimapping, and resolution

- **87 primary reads → 599 total alignments** (512 secondary): each read places at
  ~7 positions across the 6 near-identical copies. This is exactly the
  multimapping regime that naive tools mishandle (drop secondaries → miss copies;
  or count them all → inflate).
- Despite that, the reads **resolve to individual copies**. Primary placements
  distribute across the array:

| copy | primary reads | transcripts | coverage | capacity_confidence |
|---|---|---|---|---|
| c1 LOC243 | 8 | 2 | 7.9 | 1.000 |
| c2 LOC239 | 10 | 4 | 12.1 | 1.000 |
| c3 LOC240 | 14 | 3 | 10.9 | 1.000 |
| c4 LOC238 | 7 | 1 | 6.8 | 1.000 |
| c5 LOC241 | 45 | 6 | 38.0 | 1.000 |
| c6 LOC242 | 2 | 0 | 0.0 | — (below threshold) |

**5 of 6 copies are assembled.** The 6th (c6) has only 2 reads.

> ## ⚠️ CORRECTION #2 (2026-06-04) — the "collapses into ONE bundle" claim is FALSE in current code
>
> Re-measured 2026-06-04 with `RUSTLE_BUNDLE_STATS=1` on the full ampliconic
> contig (`NC_073248.2`, 11,333 reads). **The RBMY array does NOT collapse into
> one bundle.** It splits into clean per-copy bundles whose spans match the
> reference copy coordinates almost exactly:
>
> | copy | ref coords | bundle (default `-L`) | reads | tx | bundle (`--vg`) | reads | tx |
> |---|---|---|---|---|---|---|---|
> | c1 | 19,602,754–19,616,644 | 19602755–19616647 | 9 | 2 | 19599833–19616892 | 71 | 8 |
> | c2 | 19,625,715–19,639,621 | 19625716–19639637 | 9 | 4 | 19625715–19639859 | 106 | 8 |
> | c3 | 19,648,670–19,662,577 | 19648671–19662590 | 13 | 3 | 19648670–19662812 | 99 | 9 |
> | c4 | 19,671,638–19,685,525 | 19671639–19685492 | 6 | 1 | 19671638–19685780 | 97 | 11 |
> | c5 | 19,694,606–19,708,531 | 19694606–19708553 | 42 | 6 | 19694606–19731181 | 169 | 19 |
> | c6 | 19,717,578–19,730,926 | 19720059–19730114 | 2 | 0 | *(merged into c5 bundle)* | | |
>
> **Why the 2026-06-03 trace was wrong:** the secondary-bridging→mega-bundle
> mechanism only applies when secondaries are *retained*, and rustle **drops
> secondary/supplementary long-read alignments in default `-L` mode**
> (`drop_sec_supp_for_mode`, `bundle.rs:1428`). Secondaries are kept only under
> `--vg` (`vg_include_secondary`, `bundle.rs:1422`). So in the headline pipeline
> the array was *always* splitting per-copy — the collapse was a VG-mode-only
> observation, and even that no longer reproduces: in `--vg` the copies still
> mostly split (c1–c4 separate; only c5+c6 merge), and **VG family discovery
> fires** — `[VG] Discovered 15 gene family groups (78 linked bundles, 474
> multi-map reads) → 9 families after quality filter → EM reweighting complete:
> 26 reads adjusted across 6 families`. The "RBMY → 0 families" claim is stale.
>
> **Net:** rustle is NOT merging the repetitive RBMY array. The per-copy bundle
> split is faithful in both modes. (The remaining VG-side gap is over-enumeration
> from retained secondaries — c5 bundle emits 19 tx — and whether the EM families
> actually cover the RBMY bundles vs other contig families; that is a
> resolution-quality question, not a bundling-merge bug.) Genome-region
> mega-bundle scan the same day: only 2/461 chrY-amp bundles and 24/3405 chr19
> bundles swallow ≥2 gene loci, and every one matches StringTie transcript-for-
> transcript (read-through transcription both tools call, or gap-grouping where
> assembly still separates the genes). See [[project_cutoff_bundle_investigation]].

> ## ⚠️ CORRECTION (2026-06-03) — read before citing this
>
> The original framing below ("5/6 copies, each anchored `capacity_confidence
> 1.000`, the method resolves RBMY") was **wrong**, and the error is instructive.
>
> **`capacity_confidence 1.000` here is the DEFAULT, not a VG result.** It is the
> value emitted when no multimappers were apportioned — i.e. when VG copy-
> resolution **never ran**. Tracing it (2026-06-03): RBMY's 6 copies are ~23 kb
> apart, and the near-identical secondaries bridge them, so the whole array
> **collapses into ONE bundle**. VG family discovery is *inter-bundle* — it needs
> ≥2 separate bundles sharing reads — so it finds **0 families** for RBMY and the
> EM/apportionment never engages. The 5 per-copy genes were carved out by
> **ordinary position-based assembly**, not VG. So this figure did **not**
> validate copy-assignment; it showed default values. (Contrast DAZ, two distant
> *dispersed* loci → 2 bundles → 179 cross-bundle links → genuinely VG-resolved.)
>
> **What actually resolves RBMY** is the intra-bundle **tandem decomposition**
> (`RUSTLE_VG_TANDEM=1`, spec `docs/superpowers/specs/2026-06-02-intrabundle-tandem-vg-design.md`):
> it splits the collapsed mega-bundle into per-copy sub-bundles and runs the
> family graph + fingerprint-EM. On the isolated array this resolves **all 6
> copies** with *real* per-copy `capacity_confidence` reflecting each copy's
> evidence — c6 = **0.956** (the outlier, most distinguishable; recovered from 0
> transcripts), c4 = **0.220** (a core copy 99.8 %-identical to its dominant
> sibling; emitted at low confidence, honest abstain-in-place). That is genuine
> VG resolution — and the per-copy confidence *is* the scientific result.
>
> **Caveat (do not overstate):** the tandem feature is **opt-in** and currently
> fires only when the array collapses to a bundle not already swept into a
> (cross-mapping) family — i.e. on isolated/region bams. In a full-chromosome run
> RBMY's bundles get mis-assigned to an unrelated family and the pass skips them,
> so the 6-copy resolution above is from the RBMY-region bam, not yet a whole-
> genome result. The honest claim is: *VG can resolve a tandem testis array into
> per-copy isoforms with calibrated confidence; productionizing it for whole-
> genome input is in progress.*

## What this validates (corrected)

1. **The identifiability spectrum is real and the confidence tracks it.** With
   tandem mode, the 6 RBMY copies resolve at confidences that scale with their
   distinguishability: the divergent c6 (0.956) cleanly, the near-identical core
   copies lower (c4 = 0.220), and a *near-perfect* inverted duplicate (DAZ1/DAZ3
   at 99.97 %) stays non-identifiable and is flagged rather than guessed.
2. **No fabrication.** A copy is recovered only with its own observed reads;
   contested evidence yields *low* confidence, never a confident false call (the
   DAZ3 discipline).
3. **Biologically consistent.** RBMY is testis-specific; per-copy expression
   across the array is the expected result in a testis sample.

Figure: `rbmy_analysis.png` (**stale — shows the pre-correction default values;
regenerate from a `RUSTLE_VG_TANDEM=1` run**). Data: `/tmp/rbmy/`.
