# PSV↔Junction Linkage for Paralog Isoform Recovery (Layer-2 "C") — Design

**Date:** 2026-06-14
**Status:** Design approved; spec for review before planning.
**Feature flag:** `--vg-layer2-psv-linkage` / `RUSTLE_VG_LAYER2_PSV_LINKAGE` (default **off**).

## 1. Goal

Recover *harder* per-copy splice isoforms in paralog gene families by using
**within-molecule PSV→copy linkage as positive assignment evidence**: when a single
long-read molecule spans both a paralog-distinguishing position (PSV) *and* a splice
junction, the PSV tells us which copy the molecule came from, so the junction can be
assigned to that copy. Aggregating these linked assignments lets Layer-2 emit splice
forms attributed to a *specific* copy that the current per-copy assembly and Part A miss.

This is the **positive** use of allele-linkage. `decompose_family_paths` already uses
allele-linkage as a *negative* constraint (a copy may not claim a diagnostic node it has
no allele for). (C) flips it: a molecule's PSV alleles actively *assign* its junctions.

### Why this is the right lever

Paralog recovery is bounded by an **identifiability barrier**: a read is only assignable
to its copy of origin if it spans enough PSVs to beat the sequencing-error floor (~3 PSVs
/ 0.997 identity, coverage-dependent — the discriminative-position density on *exons*,
not overall %identity, is the axis). Where reads carry no discriminative signal (e.g. the
DAZ regime), nothing can resolve them and naive "secondary-at-locus = this copy's read"
logic (Part A) manufactures phantoms. (C) addresses both sides: it *uses* the PSV signal
where it exists, and *gates out* families where it doesn't.

## 2. Decisions (locked during brainstorming)

| Decision | Choice |
|---|---|
| Deliverable | **Both** — a reusable assignment engine that exposes a certificate AND drives emission of per-copy isoforms |
| FP control | **Two-level gate** — (E) family-level identifiability predictor screens whole families, then per-junction linkage (≥K linked reads, ≥N error-model PSVs) decides within identifiable families |
| Scope of v1 | **Full engine + emission + genome-wide evaluation** from the start |
| Data source for per-read PSV alleles | **Approach A** — targeted second BAM pass over family loci (reuses `bam.rs` mismatch reconstruction) |
| Relationship to Part A | **Separate opt-in channel** alongside Part A (Part A unchanged, default-on); (C) is the PSV-validated channel, opt-in |
| Default | **Off**, behind `--vg-layer2-psv-linkage`; strictly additive; VG⊇baseline floor preserved |

## 3. Architecture

(C) runs inside `run_layer2`, per family, only when `--vg-layer2-psv-linkage` is set and
`--genome-fasta` is available (PSV genotyping needs the reference). Pipeline per family:

```
family graph (copies = paths, diagnostic nodes carry per-copy alleles)
        │
   (E) identifiability gate ──fail──> SKIP family (no C emission)
        │ pass
   PSV genotyping pass (2nd BAM pass over family loci)
        │  per read: allele vector at PSV columns + junction chain
   per-read copy assignment (>=N agreeing PSVs beat error floor; else "ambiguous")
        │
   junction→copy linkage (read confidently copy X AND spans J  ⇒  vote J→X)
        │  junction J assigned to copy X iff >=K linked-read votes
   per-copy isoform assembly (assigned junctions → exon chains, chains-only >=K)
        │
   ├─> certificate (linked-read count, #PSVs, copy posterior)  → rescue report
   └─> additive emission (UnionBaseline-tagged → union-by-chain → floor-safe)
```

### 3.1 Components

**(E) `family_identifiability(fg, copies, error_rate, min_psv) -> bool`** (new, in
`vg_family/`). Computes the number of *distinguishable* PSV columns on the family's
exons from the diagnostic nodes' `per_copy_sequences` (positions where ≥2 copies differ),
and returns whether a typical read span clears the error floor. v1 rule: identifiable iff
the family has ≥ `min_psv` (default 3) distinguishable PSV columns on exonic sequence
**and** the per-family observed error rate (estimated from read mismatch rate at
non-PSV positions) leaves those PSVs distinguishable. (Refinement, documented but v1 may
approximate: weight by `P(read of median length spans ≥ min_psv PSVs)`.) Reuses
`enumerate_diagnostic_sites` (vg.rs:4225) for the diagnostic-node set.

**`PsvColumn`** (new struct): for each diagnostic exon node, the set of base positions
where copies differ, each with `per_copy_allele: Vec<(copy_id, base)>` and the genomic
position **in each copy's coordinates** (derived from the node's `per_copy_spans` + the
within-exon offset). This is the genotyping reference.

**PSV genotyping pass** (new, `genotype_family_reads(bam_path, family_loci, psv_columns, genome) -> Vec<ReadGenotype>`):
a focused second BAM pass over the family copies' genomic spans. For each read overlapping
a copy locus, for each PSV column the read covers, read the read's base at that genomic
position (CIGAR-mapped; reuses `bam.rs`'s per-position mismatch-vs-FASTA reconstruction)
and compare against the per-copy alleles. Emit `ReadGenotype { read_name_hash, psv_votes:
Vec<(copy_id, n_supporting_psvs)>, junctions: Vec<(u64,u64)> }`. The pass is bounded to
family loci × PSV columns, so it is cheap. Determinism: reads processed in BAM order;
all aggregation keyed in `DetHashMap`.

**Per-read copy assignment**: a read is **confidently assigned to copy X** iff
`psv_votes[X] >= N` (default `N=3`, env `RUSTLE_VG_LAYER2_PSV_MIN`) AND X strictly
dominates every sibling (no other copy within a margin). Otherwise the read is
**ambiguous** (not used for linkage — never guessed).

**Junction→copy linkage** (`link_junctions(read_genotypes, k) -> DetHashMap<copy_id, Vec<junction>>`):
for each confidently-assigned read, each junction it spans casts a vote J→copy. A junction
is assigned to a copy iff `votes >= K` (default `K=2`, env `RUSTLE_VG_LAYER2_PSV_MIN_LINKED`).

**Per-copy isoform assembly**: for each copy, take its linkage-assigned junctions and the
linked reads' exon coordinates, assemble exon chains (reuse `build_exons_from_chain` /
`enumerate_secondary_chains` machinery on the linked-read subset), keep chains with ≥K
linked-read support (chains-only — never a synthesized walk). Each becomes a `FamilyPath`
(reusing the existing struct) with `source = IsoformSource::PsvLinked { copy_id }` (new
variant) and a `PsvCertificate`.

**`PsvCertificate`** (new): `{ copy_id, linked_reads, n_psvs, copy_posterior }` — attached
to the emitted transcript and surfaced in the existing `RUSTLE_VG_RESCUE_REPORT` TSV.

### 3.2 How (C) differs from Part A

Part A enumerates a copy's intron chains from secondaries-at-locus **without genotyping**
— it trusts that a secondary at copy B's locus is copy B's evidence. (C) **genotypes**: a
read at copy B's locus whose PSV alleles say "copy A" is assigned to A, not B — so the
phantom Part A would emit (the M5 `u`/DAZ failure mode) is structurally prevented. (C) is
opt-in and additive; when on, it contributes PSV-validated per-copy isoforms on top of
Part A and the base recovery (union-by-chain dedups overlaps).

### 3.3 Emission & invariants

Emitted `FamilyPath`s flow through the existing path: appended to `novel_transcripts` →
`union_baseline_holdout` → union-by-chain (pipeline.rs ~19761), which is **strictly
additive** (emits a chain only if absent from VG output, exons ≥ 2; `default_min_longcov`
= 0.0 under `--vg-layer2`). A PSV-linked chain at a copy's own coordinates is a new chain →
only ADDs → **VG⊇baseline floor preserved, 0 regressions** by construction. Determinism:
`DetHashMap`/`DetHashSet`, total-order sorts, BAM-order genotyping. `--vg` (layer2 off) and
`--vg-layer2` without the new flag are **byte-identical** to today (all new work gated).

## 4. Thresholds (all env-tunable, error-model-motivated — not arbitrary)

| Name | Default | Meaning | Env |
|---|---|---|---|
| `min_psv` (N) | 3 | PSVs that must agree to confidently assign a read to a copy (the ~3-PSV/0.997 error floor) | `RUSTLE_VG_LAYER2_PSV_MIN` |
| `min_linked` (K) | 2 | confidently-assigned reads that must link a junction to a copy before it is emitted | `RUSTLE_VG_LAYER2_PSV_MIN_LINKED` |
| family `min_psv_columns` | 3 | distinguishable PSV columns required for a family to pass the (E) gate | `RUSTLE_VG_LAYER2_PSV_FAMILY_MIN` |

## 5. Edge cases / error handling

- **No genome FASTA**: (C) requires `--genome-fasta`; if absent with the flag set, log a
  warning and no-op (consistent with the rest of genome-backed Layer-2).
- **Family fails (E)**: skip — no (C) emission for that family (the phantom defense).
- **Read covers PSVs but is ambiguous** (no copy reaches N, or a tie within margin): not
  used for linkage; never assigned.
- **Read spans a junction but no PSV**: cannot be linked (the identifiability boundary) —
  contributes nothing; this is expected and silent.
- **Dispersed vs co-located copies**: (C) works per copy at the copy's own coordinates
  (genotyping is in each copy's own frame), so it does NOT depend on shared backbone nodes
  — it avoids the Part-B dispersed-paralog failure.
- **Cross-copy union span landmine**: never use `ExonClass.span` (cross-copy union); always
  `per_copy_spans` — same rule as the rest of Layer-2.

## 6. Testing

- **Unit**: `family_identifiability` (DAZ-like→fail, RABL2-like→pass); PSV genotyping
  (a read carrying copy-B alleles at copy-A's locus → assigned to B); junction linkage
  (≥K votes required); per-read ambiguity (no copy reaches N → unassigned).
- **Synthetic fixture** (`bench/sim/gen_psvlink_fixture.py`): a 2-copy family where copy A
  and copy B share an exon structure but copy A expresses an exon-skip isoform; reads of
  that isoform carry copy-A-distinguishing PSVs. The skip junction is resolvable to copy A
  **only** via the linked PSV (without linkage it is ambiguous/unassignable). Load-bearing:
  with `min_psv` honored, (C) assigns the skip to A and emits it; disabling the linkage
  check mis-assigns or drops it. Harness leg in `bench/layer2_invariant.sh`: default-off
  byte-identical; flag-on recovers the PSV-linked isoform; additive; floor holds.
- **Genome-wide** (per the chosen scope): re-run the scorecard with `--vg-layer2-psv-linkage`
  on (extend `bench/gw_threeway.sh` to add a (C)-on output column, or a parallel harness),
  and extend the attribution (`bench/gw_attribute.py`) to report **PSV-linkage net-new
  recoveries** — NCBI-corroborated transcripts that (C) adds beyond the current VG-layer
  (the 62-transcript baseline) — plus the precision delta. Standing invariant: 100%
  StringTie parity, 0 regressions, preserved.

## 7. Honest scope / limits (must remain in the spec)

- (C) fires **only** where a molecule spans both a PSV and the junction. Long HiFi reads
  make this common in identifiable families, but it is still bounded by per-read
  information; the (E) gate keeps it silent where it can't help.
- Expected contribution is **modest in absolute count** (the VG-layer's whole current
  genome-wide contribution is 62 transcripts); (C)'s value is *precision-safe* recovery of
  the harder, copy-specific splice forms in identifiable paralog families (RABL2/RBMY-like),
  with the (E) gate preventing the DAZ-regime phantoms. We measure, we don't pre-claim a
  large number.
- (C) does **not** invent novel cross-copy combinations (the homology barrier is real and
  separate) — it assigns *real, read-traversed* junctions to the correct copy.

## 8. Open questions (resolve during planning)

- Exact form of the (E) predictor: hard PSV-column count (v1) vs the
  `P(read spans ≥N PSVs)` coverage-weighted form (refinement). v1 uses the count + an
  error-rate sanity check; the refinement is a follow-up if the count proves too coarse.
- Whether the second BAM pass reuses the existing `collect_secondary_index_from_bam`
  scaffolding (a sibling collector) or a new focused reader — decided in planning by code
  fit.
- Copy-posterior computation for the certificate: simple vote-fraction (v1) vs a proper
  per-read likelihood under the error model (refinement).
