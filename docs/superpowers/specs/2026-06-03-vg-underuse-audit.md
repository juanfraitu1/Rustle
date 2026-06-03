# VG under-use audit (2026-06-03)

**Method.** A workflow fanned out 5 auditor agents across VG dimensions (env gates,
objectives, vg_hmm modules, computed outputs, SV-unique capabilities), then
**adversarially verified** each high/medium finding against the code (default to
skepticism — the burden was on proving genuine under-use, after this session's
lesson that O5 structure-borrowing *looked* unbuilt but was default-on). 28
candidates examined; ~18 confirmed genuine under-use, 3 rejected as gated-for-a-
good-reason, the rest unverified/low. **Already-addressed this session (O5 fix,
tandem decomposition, support floor, error-aware scoring, attribution benchmark)
were excluded by construction.**

An "under-use" = a VG capability that EXISTS in code but is (i) gated default-OFF
without a regression reason, (ii) computed but never consumed/surfaced, (iii)
reached only on a narrow path, or (iv) a thing VG is uniquely positioned to do but
hasn't been built out.

## Confirmed under-uses, by theme

### A. Structure / topology transfer  ⭐ top lever for the O5 "real gap"
- **`RUSTLE_VG_TOPO_BORROW`** (+ `RUSTLE_VG_TOPO_MIN_CONF`) — `transfer_assembled_topology()`
  (vg.rs:5788-5953) projects a *confident* copy's whole exon-chain onto an
  **under-assembled sister copy** via `FamilyGraph.per_copy_spans`, with safety
  rails (per-exon sister-span mapping with break-on-miss, sister-read-overlap proof
  of expression, confidence discounted ×0.5, source must be multi-exon + clear
  `min_confidence`). **Double-gated OFF; git confirms it was never disabled for a
  regression — just off-and-never-measured** (sole test is a synthetic count-non-
  decrease check; zero bench/RESULTS references). This is essentially a **built
  solution to the O5 real work item** ("each copy emitted, not partially
  assembled"). **Value HIGH, effort MEDIUM (validate, don't build).** Caveat: the
  prize size (number of expressed-but-under-assembled sisters with a qualifying
  confident source) is unquantified.

### B. Copy-number from read depth  (VG-unique, O1/O2)
- **Tandem copy-number from traversal / pileup depth** — estimate how many copies a
  fully-collapsed array holds from read depth through the repeat, vs only counting
  *assembled* copies. The substrate exists (per_copy_cov, the tandem cycle); the
  estimate is unbuilt.
- **Per-candidate-locus k-mer density** (`density_per_kb`, `n_kmer_hits` in
  positional.rs) — a copy-number proxy that is computed during the locus scan but
  not surfaced as a copy-number signal.

### C. Inversions — unify, don't just filter  (DAZ1/DAZ3)
- **Reverse-complement-aware exon-class unification for inverted paralog pairs.**
  Today inverted near-identical copies (DAZ1−/DAZ3+) are only **strand-mirror
  *filtered*** (kept apart). The FamilyGraph has no RC-aware unification, so an
  inverted pair is never modeled as one inverted-repeat structure. Building this
  would let the EM share evidence across an inversion (the DAZ case).
- **`FamilyVerdict` locus-level inversion relationship** is detected but
  under-exploited downstream.

### D. Junctions as diagnostic positions  (the O5 "structural-linkage" channel, O3)
- **Copy-specific junctions as discrete diagnostic positions** — a junction present
  in one copy and absent in a sister is a distinguishing signal that resolves reads
  SNPs alone cannot (why DAZ1−/DAZ3+ are resolvable despite near-identical
  sequence). The junction signal (`lambda_j`) exists in the EM but copy-specific
  junction *provenance* is under-used as a fingerprint axis.
- **Per-copy junction provenance on `JunctionEdge`** (which copies support each
  junction) — present in the FamilyGraph, not consumed as a read-to-copy signal.

### E. Exon-structure / TE divergence  (O3/O4, the SNP-EM blind spot)
- **`RUSTLE_VG_EXON_COV`** (+ `_ALPHA`) — `read_bundle_exon_coverage` (vg.rs:1250):
  a per-read **exon gain/loss** discriminator for copies that differ by exon count
  (4-exon vs 5-exon paralogs). Gated OFF and only reachable on the no-genome path;
  shelved as "inert" — but on the GGO panel of *structurally-equivalent* paralogs,
  i.e. a non-representative benchmark. **Never tested on its actual target.**
- **`exon_length_divergence` / copy-private (`copy_specific`) exon classes** —
  exon-length divergence and copy-private exons are computed in the FamilyGraph but
  not surfaced as **TE-derived / novel-exon** candidates (a thing VG uniquely sees).

### F. Computed-but-unused outputs
- **`abundance_min`** (EM-aware jointly-feasible lower bound on per-copy abundance) —
  computed, emitted, but not used in any downstream decision.
- **`em_weight_gap` / `em_n_sites`** (continuous per-read decisiveness) and the
  **family-level `n_decisive` / `n_moderate` / `n_uncertain`** breakdown — rich
  confidence signals computed in the EM, surfaced only in trace logs, not in the
  GTF or any gate. (These are exactly the identifiability signals the attribution
  benchmark measured — surfacing them per-transcript would make confidence legible.)

### G. Haplotype / phase
- **`--vg-phase` HP-tag haplotype bundle splitting** (`split_bundle_by_phase`) —
  partition a bundle by read HP tags; a real capability for phased multi-copy loci,
  under-used.

### H. Novel-copy rescue (partial — read carefully)
- **k-mer-primed novel-copy rescue** (`run_rescue_with_bundles` → `synthesize_bundles`)
  and **profile-HMM scoring against candidate loci ("ghost copy" injection)** are
  genuine capabilities, but **the full profile-HMM rescue path (`fit_profiles`) and
  `vg_scan_novel_loci` were verified as gated-for-a-good-reason** (HMM-EM is RETIRED;
  `do_hmm` hard-coded false). So treat the HMM-dependent parts as *correctly* off; the
  k-mer/positional parts are the live, under-used portion.

## Rejected (verified gated-for-a-good-reason — do NOT re-enable blindly)
- **`RUSTLE_VG_FAMILY_RESCUE` / `RUSTLE_VG_FAMILY_TOLERANT`** — cross-bundle family
  chain-rescue; off for a reason.
- **`vg_scan_novel_loci`** — genome-wide positional novel-paralog scan; off for a reason.
- **Full profile-HMM novel-copy rescue (`fit_profiles`)** — HMM-EM retired.

## Dangerous-but-real (handle with the DAZ3 discipline)
- **`RUSTLE_VG_FAMILY_BOOST`** (vg.rs:2029) — the ONLY mechanism that recovers
  **pure-secondary copies** (primary≈0 — *exactly DAZ3's echo pattern*). Genuine
  under-use for O2, NOT superseded by the session's TandemCopy-scoped support floor.
  **But:** its "+7 transcripts" headline is a raw count delta, never precision-
  validated, and boosting a primary=0 copy is the fabrication risk that makes DAZ3 a
  false positive. **Do not enable without first running the never-run precision check.**

## Recommended order  (REVISED after the TOPO_BORROW validation below)
1. ~~Validate `RUSTLE_VG_TOPO_BORROW`~~ → **DONE; it is INERT** (see validation section).
   The real lever it exposed: **divergence-robust cross-copy exon unification** in the
   FamilyGraph (`merge_singletons_by_sequence`, k-mer Jaccard collapses at 95–99%
   identity → single-copy ExonClasses → all structure-sharing dead). Fixing THIS
   unlocks TOPO_BORROW + completion + borrow together. **New #1.**
2. **Surface the confidence outputs** (theme F: `em_weight_gap` / decisive breakdown
   per-transcript) — cheap, makes identifiability legible, complements the benchmark.
3. **Copy-number-from-depth** (theme B) — VG-unique, modest build, strengthens O1/O2.
4. **FAMILY_BOOST precision check** — high value, but gated on the DAZ3 fabrication
   risk; never enable without the precision measurement.

## VALIDATION: RUSTLE_VG_TOPO_BORROW is built but INERT (2026-06-03)

Ran the audit's #1 recommendation — validate TOPO_BORROW on a real paralog family.
Result: **it fires on nothing.** Tested three ways:
- chrY `ggo_Y.bam` 28 families (`--vg`): 316=316 transcripts, **0 topology_borrow**, no `[VG-TOPO]`.
- RBMY tandem family: **0 topology_borrow**.
- Controlled synthetic — copy 0 fully assembled (`copy_confidence 0.739` ≥ 0.7 → a
  qualifying source), copy 2 deliberately starved (3 reads → no transcript, the
  ideal under-assembled sister): **still 0.**

**Definitive root cause** (added `RUSTLE_VG_TOPO_TRACE`, commit 8a5d6e6): the source
qualifies (`n_src_txs=1, n_copies=3`) but projecting its 5 exons onto each sister
gives `projected=0` — **the first source exon maps to no ExonClass shared with the
sister.** Confirmed via COMPLETION_TRACE: **every ExonClass has `max_sib=0.00`** —
no sibling copy contributes to any exon class. So at 97% identity
`merge_singletons_by_sequence` (build_family_graph `merge_min_jaccard=0.30`,
vg.rs:480) **does NOT unify the dispersed copies' homologous exons into shared
ExonClasses** — `per_copy_spans` stay single-copy. k-mer (k=15) Jaccard collapses
fast with SNPs (a 3%-divergent 300 bp exon has ~9 SNPs disrupting ~135 of ~285
15-mers → Jaccard ≈ the 0.30 boundary), so 95–99% paralog pairs fall below the bar.

**This single bottleneck gates the ENTIRE structure-sharing theme.** TOPO_BORROW,
the completion-node injection, and coverage/junction borrowing ALL project via
shared `per_copy_spans` / `max_sib`. When the exon-unification fails (moderate
divergence), all three are inert. It works for RBMY (99.8% → exons merge,
`max_sib>0`, completion engages — as this session's tandem fix showed) but fails
for the exact O5 "real gap" case (expressed paralog pairs, 95–99%, dispersed).

**REVISED recommendation.** TOPO_BORROW is NOT the easy validate-and-ship win — it
is *correct code aimed at a target the FamilyGraph never delivers*. The genuine,
high-leverage O5 work is **divergence-robust cross-copy exon unification**: replace
or augment the k-mer-Jaccard singleton merge with an alignment/minimizer-based
homology merge (or lower/auto-scale `merge_min_jaccard` with a spurious-merge
guard), so homologous exons across 95–99% copies share an ExonClass. That single
fix unlocks TOPO_BORROW + completion + borrow simultaneously for the real O5
paralog-pair case. (This is the same shape as the session-opening lesson: a
capability that *looks* like the answer but doesn't engage — validate that it
FIRES, not just that the code is correct.)
