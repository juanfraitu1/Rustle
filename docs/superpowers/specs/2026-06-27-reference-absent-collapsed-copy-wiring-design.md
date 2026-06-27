# Reference-Absent (Collapsed) Copy Wiring — Design Spec

**Date:** 2026-06-27
**Objective:** O4 — find copies that exist only in the reads, not the genome assembly.
**Scope of THIS spec:** the **collapsed** reference-absent route only (≥2 paralogs the assembly/minimap2
merged onto one reference locus). The **divergent** route (paralogs that map nowhere → de-novo assemble
unmapped reads) is explicitly **deferred** to a later milestone.

> Status: design for review. Derived from the 2026-06-27 multi-agent design analysis (workflow
> `wf_75e6c205-c0c`) and verified against the current code. Not yet a plan; do not implement until approved.

---

## 1. Goal

Wire the **already-built read-covariation copy discovery** (`split_locus_copies` →
`CopyIsoform{allele_vector}`) into the **assignment path** so a collapsed paralog that is *absent from the
assembly* becomes a first-class copy reads can be assigned to — **without regressing the frozen O2 headline**
(75.1% assigned / 24.8% certified-tied / 99.9% of decisive / silver 99.9%), and **without fabricating a copy
out of a diploid allele, RNA edit, or alignment artifact.**

The unifying idea (sound, and the target architecture): a reference-PRESENT copy and a reference-ABSENT copy
are *the same object* — a path through one family variation graph whose bubbles are PSV columns. The **only**
difference is how the copy's `.seq` is obtained: reference-fetched (`build_spliced_seq`) vs read-consensus.
Once an absent copy is a synthetic `DenovoTranscript` appended to the copy slice, `build_family_profiles` /
`discover_psvs` / `assign_read_editing` treat it identically. **Empty absent-set ⇒ byte-identical O2.**

## 2. Background — why the two PSV sources are disjoint today (and why that makes O2 safe)

- **Ref-present PSVs** (`discover_psvs`, `copy_assign_pipeline.rs:234`): align each copy's spliced **reference**
  sequence to `copies[0]` (poasta). Non-circular: the variant calls never touch a read assignment.
- **Read-covariation discovery** (`split_locus_copies` / `split_readchain_by_psv`, `copy_split.rs:334/45`):
  partition a locus's reads into mutually-exclusive co-varying PSV haplotypes → `CopyIsoform { intron_chain,
  allele_vector: Vec<Option<u8>>, read_count, identifiable }` (`copy_split.rs:28`). Non-circular for the
  **existence** claim: co-variation is a property of the read multiset, not of any assignment.
- **The gap:** `recover_collapsed_copies` (`denovo_pipeline.rs:367`) calls `split_locus_copies` but returns a
  **`usize` count** and **discards** the `CopyIsoform`s. `build_family_profiles` (`copy_assign_pipeline.rs:324`,
  takes `&[&DenovoTranscript]`) is the **sole** production `CopyProfile` source and never consumes split
  output. So the two sources are fully disjoint → today an absent copy can never enter assignment → O2 is
  trivially safe. This spec closes that gap **under guards strong enough to keep O2 safe by construction.**

## 3. The hard boundary that shapes the output (non-negotiable)

Observed expressed-haplotype multiplicity = f(n_loci × ploidy × allele-specific-expression) — a **many-to-one,
non-invertible** map. From RNA reads alone you can prove *"≥2 distinct expressed entities here"* but **not**
*"this entity is a COPY rather than a diploid allele / ASE haplotype."* Depth doesn't disambiguate (depth = CN
× expression × ASE). This is an **information boundary requiring DNA / parental copy-number**, corroborated by
the codebase's own results (hidden-copy real headroom ≈0 on GGO; O4 FP bound 7.39% het-dominated; the 4
"confirmed" absent copies are all MHC — the worst het regime). The K=0 frontier (copies identical over exons
*and* splice sites) is unconditionally unresolvable from RNA — detect/flag/co-quantify only, never per-read.

**Consequence — two output states, never one forced "copy":**
1. **COPY** — admitted into assignment, only when ALL admission discriminators (§4.4) pass.
2. **DNA-needs / parCN candidate** — emitted, flagged, co-quantified, but **NOT** added to the copy set (does
   not affect assignment). This is a first-class output, not a footnote.

## 4. Design

### 4.1 Two-stage freeze (the O2-protection core — mandatory, do first)

Single-stage addition of an absent copy **can regress O2**: raising `n` tightens the Bonferroni threshold
`thr = alpha/(n-1)` (`copy_assign.rs:338`) **and** adds a competitor, so `min_p` and `p_read` can only rise →
a correctly-Assigned read can demote to Tied/Ambiguous (**recall loss**); and a phantom copy can become argmax
for a uniquely-mapped read truly from a reference copy, flipping `best_copy != mapped_copy` (**silver loss**).
So additivity is **precision-safe-per-read but recall-fragile and silver-unsafe** unless we freeze:

- **Stage 1 — reference-only assignment.** Run `assign_family_detailed` exactly as today (ref copies only).
  This reproduces O2 **byte-identically** (it IS today's path). Record each read's Stage-1 status.
- **Stage 2 — absent-aware re-threading, restricted.** Build the augmented copy set (ref + admitted absent
  copies) and re-assign **only**: (a) reads Stage-1-classified **Tied or Ambiguous** (the abstain pool), and
  (b) reads at **single-ref-copy collapsed loci** (where Stage 1 had no PSVs at all). **Freeze** every
  Stage-1 **Assigned** read at a multi-ref-copy family — its result is carried through unchanged.
- **Guarantee:** the O2 *assignment set* is a literal subset of the output by construction. Recall can only
  grow (abstained reads may now resolve). **But aggregate silver and %-assigned must be MEASURED** (§6) — the
  freeze does not by itself prove silver non-decreasing, because Stage 2 can still misassign a newly-captured
  unique-mapper at a single-ref-copy locus to a phantom copy.

### 4.2 The consensus-`.seq` + synthetic-`DenovoTranscript` primitive (new code)

`CopyIsoform` has `intron_chain` + `allele_vector`, **no `.seq`**, and `consensus_haplotype`
(`copy_split.rs:394`) returns only the allele vector — none of these is appendable to `&[&DenovoTranscript]`
(which `discover_psvs` requires to satisfy the `exon_maps[c].len() == seq.len()` invariant). New primitive:

`fn collapsed_copy_to_transcript(iso: &CopyIsoform, locus_ref: &DenovoTranscript, genome) -> Option<DenovoTranscript>`

- **Recommended (v1, substitution-only):** take the locus's spliced **reference** sequence (via
  `iso.intron_chain` + `build_spliced_seq`) and **overlay** `iso.allele_vector` at the PSV column positions →
  a synthetic spliced seq carrying the copy's distinguishing bases, with an `exon_map` derived from the same
  intron chain (so the `len == seq.len()` invariant holds by construction). Cheap, exact for substitutions,
  blind to copy-private indels.
- **Alternative (v2, full):** POA consensus over the copy's read group → spliced consensus (captures indels).
  Heavier; defer unless v1's indel-blindness proves limiting.
- Unit-test the invariant `transcript.seq.len() == exon_map(&transcript).len()` and that an empty/degenerate
  allele vector yields `None` (never a malformed copy).

### 4.3 Cross-frame seeding (required for the admission test)

`select_identifiable_copies` de-dups candidates **within the pileup frame** only. Reference copies' allele
vectors live in the `copies[0]`-offset frame. To ask "is this read-cluster ≥K-distinct from every *reference*
copy?" (and to catch a contaminating paralog B leaking into locus A), the reference copies' alleles must be
brought into the pileup frame via the existing ref-vs-ref alignment. This cross-frame seed is **new** (the
prior design wrongly assumed it was ~present). It is the bridge that makes §4.4's distinctness test meaningful.

### 4.4 Admission gate — a `CopyIsoform` becomes a COPY only if ALL hold

The per-read significance gate guards **precision but NOT copy-count inflation** (a het-allele admitted as a
copy confidently splits one real copy into two, each with confidently-assigned reads — the per-read bound is
blind to this). Copy-count inflation is controlled **only** at discovery, by requiring ALL of:

1. **≥K-distinct from every reference copy** in the shared frame (§4.3) — else it's a known ref copy, not absent.
2. **Multiplicity / structure:** ≥3 mutually-exclusive co-varying clusters at the locus **OR** a co-segregating
   **copy-private splice junction**. (A single biallelic diploid het yields at most 2 clusters; ≥3 rules it
   out. A private junction is positive copy evidence alleles can't fake.)
3. **Strand-symmetric substitution spectrum** in transcription orientation — RNA editing is A→G (T→C on minus)
   only; a true copy's divergence is symmetric. Port the editing test into **discovery** (`discover_locus_psvs`),
   not only the assignment-path `detect_editing_columns` (`copy_assign_pipeline.rs`).
4. **Genome-remap gate:** the copy consensus maps back to a reference locus at **<98% identity** within a
   **3–20% divergence band** (the `promote_hidden_copies.py` test, brought in-pipeline). This is **non-circular**
   (a fresh falsifiable query) and **defuses the MAPQ-0 leak**: a leaked paralog-B read consensus maps ~99.7%
   to locus B and fails the <98% gate.
5. **Not editing-flagged** and **passes an indel-proximity / low-complexity / base-QV mask** before clustering
   (the DSFAM213 correlated-alignment-error failure mode).

**Fail any ⇒ DNA-needs/parCN output state (§3), not a copy.** A 2-cluster substitution-only locus is
copy-vs-allele unresolvable from RNA → always DNA-needs.

### 4.5 Anti-circularity for the certificate (the subtle one)

Building the absent copy's consensus **from** the cluster reads and then scoring **those same reads** against
it makes the IsoCon certificate **anti-conservative for the discovered copy's own reads** (the `alpha`
guarantee assumes copy profiles are independent of the reads — true for reference-fetched `.seq`, false for
read-consensus `.seq`). Mitigation (pick one, v1 = the simplest that's honest):
- **Hold-out (preferred):** build the consensus on one half of the cluster, certify on the disjoint remainder.
- **Flag:** mark assignments to a discovery-coupled copy as `discovery_coupled=true` so downstream
  metrics/quant can treat them separately. (At minimum, v1 MUST flag; hold-out can follow.)

## 5. Out of scope (this milestone)

- **Divergent route** (assemble unmapped/poorly-mapped reads → contig → realign reads to contig). Reads carry
  genome-coordinate CIGARs; a contig-coordinate copy needs a new feature-extraction path
  (`fill_psv_obs`/`read_features` read in the copy's frame). Separate milestone; do not block collapsed
  delivery on it.
- **Layer-2 `--vg` indel composition** (`psv_columns_for_family` handles equal-length per-copy seqs only).
- Changing the ref-present PSV path, the significance gate math, or the editing filter semantics.

## 6. Acceptance — O2 non-regression measurement protocol (the real gate)

The two-stage freeze guarantees the assignment *set* is an O2 superset, but the **headline metrics must be
measured, not assumed**:

1. Run the co-located O2 substrate (the 74-family `o2_regions.txt`) **with absent-copy wiring OFF** → must be
   byte-identical to the current `o2_definitive.assignments.tsv` (regression guard on the Stage-1 path).
2. Run **with wiring ON** and require: **silver non-decreasing**; **%-assigned non-decreasing**; the assigned
   set a near-superset of O2 (report every O2-Assigned→Tied/Ambiguous flip — should be ∅ given the freeze).
3. Report: # absent copies admitted vs # routed to DNA-needs; # reads newly resolved from the abstain pool;
   any `best_copy != mapped_copy` silver flips at single-ref-copy loci (the freeze's blind spot).
4. Sim5x labeled-truth check: inject a known collapsed copy; require it is admitted and its reads correctly
   assigned, and that a planted **diploid het** (not a copy) is routed to DNA-needs, not admitted.

Ship only if (1) is byte-identical and (2)/(4) hold.

## 7. Components / files

| Component | File:line | Change |
|-----------|-----------|--------|
| Collapsed discovery | `copy_split.rs:334` `split_locus_copies` → `CopyIsoform` | reuse; add discovery discriminators (§4.4.2/3/5) into `discover_locus_psvs` |
| Discard→emit | `denovo_pipeline.rs:367` `recover_collapsed_copies` | return the `CopyIsoform`s (not a `usize`) |
| New primitive | `copy_split.rs` (new) | `collapsed_copy_to_transcript` (§4.2) |
| Cross-frame seed | `copy_split.rs` `select_identifiable_copies` | seed with reference allele vectors mapped into the pileup frame (§4.3) |
| Admission gate | new (discovery side) | the all-of gate (§4.4) + DNA-needs output state |
| Remap gate in-pipeline | port from `bench/promote_hidden_copies.py` | §4.4.4 |
| Two-stage freeze | `denovo_pipeline.rs:~507` (the colocated-family assign loop) + `assign_family_detailed` | §4.1 |
| Certificate hold-out/flag | `copy_assign.rs` / `copy_assign_pipeline.rs` | §4.5 |
| Existing guards reused | `hidden_copy.rs:76`, `detect_editing_columns` | §4.4 |

## 8. Decisions (locked for v1, 2026-06-27)

- **consensus-`.seq` = reference+allele OVERLAY** (§4.2 v1): spliced reference at the locus with
  `iso.allele_vector` overlaid at the PSV columns. Substitution-only; a collapsed copy suspected of a private
  indel is **flagged** (`indel_suspected`) and routed to DNA-needs rather than mis-built. Full-POA consensus
  (v2) deferred.
- **K (distinctness) = `min_p`-derived**, not a fixed column count: a candidate is "≥K-distinct from a
  reference copy" iff, treated as a read, its identifiability bound `min_p` against that copy is `< α` (the
  same threshold-light, Canzar-aligned certificate the assignment gate uses). The `hidden_copy` ≥12-column bar
  is **not** used for admission.
- **Certificate = FLAG in v1** (§4.5): assignments to a discovery-coupled copy carry `discovery_coupled=true`;
  the certificate is computed as usual but disclosed as anti-conservative for that copy's own reads. Hold-out
  (build-on-half / certify-on-remainder) is a v2 follow-up.
- **DNA-needs output = a new `<out>.dna_needs.tsv`** (one row per rejected candidate: locus, cluster count,
  why-rejected, divergence, read support), separate from the assignments/quant outputs.

Remaining genuinely-open (decide during implementation, not blocking): whether the cross-frame seed (§4.3)
reuses `discover_psvs`'s existing ref-vs-ref alignment output directly or recomputes a lightweight map.

## 9. Risks (and the guard that catches each)

| Risk | Guarded by | Residual |
|------|-----------|----------|
| Diploid allele admitted as a copy (dominant FP) | §4.4 discriminators + DNA-needs state | 2-cluster sub-only → DNA-needs (correct, not a copy) |
| Certificate anti-conservative on own reads | §4.5 hold-out/flag | flag-only v1 still anti-conservative; disclose |
| RNA editing → fake copy-private column | §4.4.3 strand-symmetry in discovery | A→I only; other editing types uncovered |
| Correlated alignment error (homopolymer/STR) | §4.4.5 indel/QV mask | needs the mask actually ported |
| O2 recall/silver regression | §4.1 freeze + §6 measurement | single-ref-copy locus silver flip (measured, §6.3) |

---

### Self-review
- Placeholders: none (all "TODO"-shaped items are explicit Open Questions §8, scoped Out-of-Scope §5).
- Consistency: the two output states (§3) thread through §4.4 and §6.4; the freeze (§4.1) and the measurement
  (§6) are the same O2-protection contract stated from design and acceptance sides.
- Scope: collapsed route only; divergent + VG-indel explicitly deferred (§5).
- Ambiguity: the consensus-`.seq` method and K threshold are the two genuine forks — captured as Open
  Questions with a recommended v1 default, not left implicit.
