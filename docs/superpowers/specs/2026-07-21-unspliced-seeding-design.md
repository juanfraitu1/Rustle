# Position-Aware Seeding of Unspliced Reads — Design

**Date:** 2026-07-21
**Status:** approved design, pre-implementation
**Motivation.** The advisor's over-merge objection traces (via two rounds of validation) to a
**detection-seeding bug**, not a merge bug. The 28 "distinguishable-but-merged" Soto members —
pseudogene/retrocopy loci whose reads are largely **unspliced** — lose their reads before any
copy decision is made.

**Root cause (verified in code).** `pass1_skeletons_robust` (`denovo_assemble.rs`) groups reads
into `BTreeMap<(chrom, introns), …>`. Unspliced reads have `introns = []`, so **every unspliced
read on a chromosome lands in the single key `(chrom, [])`**. The resulting skeleton's `start` is
the chromosome-wide k-th-smallest start and `end` the k-th-largest end — a 300–680 kb span that
exceeds `MAX_SPLICED = 300_000` and is rejected by `assemble_gate`. The member's read evidence is
deleted at seeding. The codebase **already documents this exact hazard** at
`denovo_assemble.rs:685` ("keying on the empty intron chain would union every unspliced read on a
chromosome") — but the fix was applied only to `locus_support`/`thin_loci`, never to
`pass1_skeletons_robust` itself.

**Not** `distinct_locus_reps` (the prior effort's target — requires literal overlap), **not**
`collapse_loci_span_aware`, **not** O2 assignment. Confirmed by tracing 3 cases end-to-end.

---

## 1. The change (one function)

In `pass1_skeletons_robust`, split the read stream by splice status:

- **Spliced reads** (`introns` non-empty) → the current `(chrom, introns)` keying, **untouched**.
  Byte-identical for any spliced-only input.
- **Unspliced reads** (`introns` empty) → per chromosome, **single-linkage overlap clustering**
  on `[ref_start, ref_end]`: reads whose spans overlap join one cluster; each cluster becomes its
  own `Skeleton { introns: [], start: cluster_min_start, end: cluster_max_end, n_reads }`,
  kept only if `n_reads >= min_reads`. The k-smallest-start / k-largest-end robust-boundary logic
  is applied **per cluster**, not per chromosome.

This stops the whole-chromosome pooling: each unspliced member seeds its own skeleton within its
own extent (well under `MAX_SPLICED`), so its reads survive to the copy-decision stages.

### Why single-linkage overlap (threshold-free)

Overlap needs no distance constant. Long HiFi reads tiling one ~1–3 kb pseudogene overlap heavily
→ one cluster. Copies 60 kb–2.5 Mb apart never overlap → separate clusters. This is the same
single-linkage-by-span merge `thin_loci` (`rescue_pipeline.rs:67`) already uses — no arbitrary
number enters. **Pure overlap, no gap tolerance** (long reads rarely leave gaps within a short
pseudogene; a gap tolerance would introduce exactly the arbitrary threshold this avoids).

---

## 2. What happens downstream (identifiability still governs)

The fix only **restores the reads**; it does not decide copy identity. Once an unspliced member
seeds its own skeleton, the existing machinery takes over:
- **Distinguishable members** (the 28, with 11–40 unique-mapper reads) → recovered as separate
  copies (member sensitivity rises).
- **True-K=0 members** (the 36 "exon-homogenized") → seeded as loci, but the read-conflict /
  significance machinery **abstains** on assignment (χ(H) ties them) — the honest "copy exists,
  unresolvable from RNA" state, not a fabricated distinction.

The fix and the identifiability theorem compose: seeding gives every member a fair chance to be
seen; χ(H) still decides what can be resolved.

---

## 3. Components & boundaries

| Unit | Responsibility | Input | Output |
|---|---|---|---|
| `cluster_unspliced(reads, min_reads, k)` | single-linkage overlap clustering of empty-chain reads → skeletons | the unspliced reads | `Vec<Skeleton>` (introns=[]) |
| `pass1_skeletons_robust` | seed skeletons; spliced path unchanged, unspliced path delegated to the clusterer | all reads | `Vec<Skeleton>` |

`cluster_unspliced` is a pure, independently-testable function (sort by start; single-linkage
merge overlapping spans; per-cluster robust boundaries + `min_reads` filter). The spliced branch
of `pass1_skeletons_robust` is byte-for-byte the current code.

---

## 4. Validation (result-changing — old-vs-new on Soto)

Same old-vs-new per-region method the prior effort used, now expected to actually move the 28.
Success criteria:

1. **The 28 seed/recover.** For the highest-unique-read merged cases (ID_26/40, ID_261/24,
   ID_260/11 …), the NEW binary forms a skeleton at the member's locus where the OLD binary
   deleted it; member sensitivity rises above **276/362 (76.2%)**.
2. **Precision held at 100%.** Every newly-seeded locus overlaps a real annotated Soto member —
   no spurious noise splits (the `min_reads ≥ 3` floor is the guard). If any new locus does NOT
   correspond to an annotated member, that is a precision regression → stop and investigate.
3. **The 36 K=0 behave honestly** — seeded-but-tied or still merged, never a *false* extra copy.
4. **Spliced-only regions byte-identical** (spliced path untouched) — spot-check a spliced Soto
   family's `families.tsv` md5 old-vs-new.
5. **`cargo test` green**, plus unit tests: `cluster_unspliced` separates two non-overlapping
   unspliced read groups into two skeletons and merges one overlapping group into one; a
   whole-chromosome scatter of unspliced reads at 3 distinct loci yields 3 skeletons, not 1
   giant rejected one.
6. **The 10 "unseeded"** are re-checked (`soto_floor_decompose.py`) — likely improved (same root
   cause), reported as a bonus, not required.

---

## 5. Risks

- **Noise clusters → spurious loci.** A handful of stray unspliced reads could form a cluster.
  Mitigation: the existing `min_reads` (≥3) floor per cluster; the precision check (new loci must
  map to annotated members) is the backstop.
- **Single-linkage daisy-chaining.** If unspliced reads formed a continuous overlapping chain
  across a large region, single-linkage could over-merge. Discrete pseudogene copies do not chain
  (they are separated by non-transcribed genomic distance); validation confirms the 60 kb–2.5 Mb
  gaps produce separate clusters.
- **Downstream count shifts on K=0.** Newly-seeded K=0 members change what the catalog reports
  (loci now present, assignment abstains). This is correct behavior; the precision check ensures
  the new loci are real members, not fabrications.
- **Performance.** Per-chromosome single-linkage over unspliced reads is O(n log n) (sort +
  linear merge) — negligible vs the existing pipeline.

---

## 6. Success criteria

- `pass1_skeletons_robust` seeds unspliced reads position-aware; spliced path unchanged.
- `cluster_unspliced` unit-tested (separate/merge/multi-locus-scatter).
- Soto member sensitivity **> 76.2%**, precision **= 100%**, the 36 K=0 members produce no false
  copies, spliced-only regions byte-identical.
- The 28 merged cases traced in diagnosis now seed at their loci (old-vs-new).
- `cargo test` green; scope limited to the one function + its clusterer.
