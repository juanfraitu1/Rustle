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

## Exon-unification fix attempt (2026-06-03, commit 749a607)

Pursued the revised #1 lever. The chase went deeper than "lower the merge threshold":
- The minimizer-Jaccard bar (`merge_min_jaccard=0.30`) IS too strict at 95-99%
  identity, but lowering it (env `RUSTLE_VG_FAMILY_MERGE_JACCARD`, added) only lifted
  shared ExonClasses 0→6/15 and hit a ceiling — because most exons weren't even
  reaching the merge.
- **Real bug found:** `merge_singletons_by_sequence` split clusters by EXON COUNT
  (`len>1` ⇒ "multi-copy, keep as-is"), but for DISPERSED copies a `len>1` cluster is
  usually **same-copy isoform variants** grouped by position — a single-copy cluster
  wrongly treated as already-unified, so its homologous core exon never merged across
  copies. Fix (`RUSTLE_VG_FAMILY_MERGE_SAME_COPY`, opt-in): route clusters by distinct
  **copy** count, not exon count. Result on the starved synthetic: **shared
  ExonClasses 6/15 → 15/27** — cross-copy unification more than doubled.
- **Still multi-layer:** even with 15/27 shared, TOPO_BORROW emits 0 — a FURTHER
  projection blocker remains (the assembled-source-exon ↔ family-graph-span match at
  vg.rs:5852, likely SPAN_TOL / boundary mismatch, before the cross-copy check). So
  recovering a starved copy end-to-end is genuinely 3 layers: (1) merge threshold,
  (2) same-copy-cluster routing [both addressed, opt-in], (3) source-exon↔span match
  [open].

**Status:** the merge fix is a real improvement to FamilyGraph completeness
(benefits completion/borrow broadly), committed opt-in (default byte-identical: DAZ
unchanged, RBMY c4/c6 retained, de-novo untouched). The full TOPO_BORROW-fires goal
needs layer (3) — a well-scoped next step (instrument the per-exon src-side match in
`transfer_assembled_topology` to see whether boundaries miss within SPAN_TOL).

## Layer 3 RESOLVED + Layer 4 revealed: TOPO_BORROW projection fixed, but emits 0 on real data BY DESIGN (2026-06-03)

Continued the layer-3 debug. Built a controlled starved synthetic (3 dispersed copies
200 kb apart at 97% exonic identity; copy-0 covered with 40 reads, sisters 1/2 starved
to 2 reads each, linked by copy-0's secondary alignments) and instrumented the per-exon
projection in `transfer_assembled_topology` to split the failure modes.

**Two real projection bugs found and fixed:**
1. **`find` → `find_map` (correctness fix).** The projection used `fg.nodes.iter().find()`
   to locate the ExonClass containing the source exon, then checked THAT node for a sister
   span. But a source exon can belong to BOTH a copy-private node and a shared node; `find`
   stopped at the first (often the private one) and failed (the old "MODE_B"). `find_map`
   searches on for a node shared with the sister. On the synthetic this alone projected all
   5 exons. **Only active under `RUSTLE_VG_TOPO_BORROW` (opt-in); default path untouched.**
2. **Offset bootstrap (`RUSTLE_VG_TOPO_OFFSET`, opt-in).** For exons with genuinely no
   shared ExonClass (k-mer Jaccard ≈0 for short, divergent exons — threshold-insensitive:
   merge at 0.05/0.02/0.0 all fail identically), learn the colinear copy-to-copy shift from
   the exons that DID unify (median of per-class src↔sister start deltas) and project the
   rest by that shift. On real GOLGA8 this projected **4287 exons** that the class lookup
   could not — a genuine capability gain over the prior totally-inert state.

**Verified end-to-end on the synthetic:** with the projection fix + a diagnostic
`RUSTLE_VG_TOPO_FORCE_EMIT` (nominal coverage + phantom-guard exemption), borrowing recovers
both starved sisters at the exact truth coordinates → **gffcompare 100% Sn / 100% Pr, 3/3
loci**. So the projection machinery is now correct.

**Layer 4 (the real wall): TOPO_BORROW still emits 0 transcripts on real data — and that is
CORRECT, not a bug.** Two safety layers, both verified firing on GOLGA8/chrY:
- **`has_reads` (sister-read-coverage check, upstream of emission).** Of GOLGA8 fam=11's
  projections, the 3 copy-pairs with a learnable offset projected fully, then ALL failed
  `has_reads`: the borrowed structure lands where the real sister has no reads, so it
  refuses to fabricate. `FORCE_EMIT` does NOT resurrect these (it only bypasses the later
  phantom guard) — GOLGA8 + FORCE_EMIT still emits 0. chrY/DAZ: same (projects 8 exons via
  offset, then `has_reads`-abstains).
- **Copy-support phantom guard (post-emission).** A borrowed copy has, by construction,
  ~0 independent coverage / sub-τ independent support — the exact phantom signature the
  guard suppresses (DAZ3). On the synthetic this is what removed the +2 until FORCE_EMIT
  exempted them.
- Also structural: **9 of 12 GOLGA8 copy-pairs share zero ExonClasses**, so the offset
  can't even bootstrap (585 exon SKIPs).

**Conclusion (honest, thesis-relevant).** Recovering a starved copy end-to-end only succeeds
when the starved copy has genuine read coverage over the borrowed structure (the synthetic).
A truly ~0-read copy can only be emitted by fabrication — forcing past BOTH the `has_reads`
check and the phantom guard — which is precisely the DAZ3 false-positive mechanism. The O5
"real gap" (≈0-coverage copies) is therefore structurally non-recoverable without fabrication;
the identifiability boundary is real, not a missing feature. This vindicates the multi-layer
safety design. TOPO_BORROW's projection is now correct and fires on real data; its (correct)
silence at output is the anti-fabrication guards doing their job.

**Shipped (all opt-in, default byte-identical — every change is inside the
`RUSTLE_VG_TOPO_BORROW` path or its own env gate):** find_map projection fix,
`RUSTLE_VG_TOPO_OFFSET` bootstrap, `RUSTLE_VG_TOPO_FORCE_EMIT` diagnostic (clearly labeled
as the fabrication-measurement path, NOT a default), MODE/via trace refinement.

### Correction (same session): the offset is synthetic-exact but real-data-IMPRECISE — do not overclaim

Verified the real-data projections by contig (the earlier check used the wrong contig).
On chrY the source is RSTL.13 / FAM_34 with `family_verdict "not_expressed_here"`,
`family_n_expressed "0"` — i.e. the sister copy is **genuinely unexpressed**, nothing to
borrow. The +9.8 Mb median offset projected the sister to ~26.07 M, but the nearest real
genes on that contig sit at 26.37 M+ — the single median shift landed **~300 kb off**.
So:
- The "4287 exons projected on GOLGA8" figure is **inflated and not a capability win**:
  it is ~163 redundant source transcripts from one copy × their exons, all projected by
  one crude median offset to coordinates the sister has no reads for, then (correctly)
  dropped by `has_reads`. A single affine offset cannot model real paralog
  indels/rearrangements; it is exact only on the colinear synthetic.
- The honest, durable results are: (1) the **`find_map` correctness fix** (real bug,
  verified on the synthetic, default-safe), and (2) the **complete characterization of
  why TOPO_BORROW emits 0 on real data** — the candidate sisters are genuinely
  unexpressed, and the `has_reads` + phantom-guard layers correctly refuse to fabricate.
- The offset bootstrap is therefore a **synthetic-only demonstrator**, not a real-data
  lever. It stays opt-in and clearly labeled; it should NOT be presented as recovering
  real copies.

Net: TOPO_BORROW's projection path is now correct and *can* fire, but there is **no real
GOLGA8/chrY beneficiary** because the "starved" sisters are genuinely not expressed. The
O5 ~0-coverage gap remains structurally non-recoverable without fabrication. No fabricated
win; the offset's real-data value is honestly ~nil.

## Borrow-beneficiary hunt: NONE found in chrY or GOLGA8-region (2026-06-03)

Asked directly: is there a real paralog family where TOPO_BORROW would *correctly* complete
an under-assembled-but-genuinely-full-structure copy (a legit beneficiary, not fabrication)?
Scanned both available multi-copy datasets. **No genuine beneficiary exists**, for a
structural reason, in two regimes:

- **chrY (ggo_Y.bam):** every discovered family assembles exactly ONE copy; 61/62 VG
  transcripts report `family_n_expressed 0`. The sister copies are genuinely unexpressed
  (~0 reads). Nothing to borrow.
- **GOLGA8-region (golga8_region.bam):** families assemble MANY copies, often
  over-enumerated (FAM_11 copy 3 = 163 isoforms). The only "thin" copies (FAM_1 copy 0:
  4 exons, copy 3: 3 exons vs sibling copy 1's 24) are genuinely divergent paralogs, NOT
  under-assembled colinear copies:
  - copy 3's locus spans **250 kb with 3 exons** while rich sibling copy 1 is **70 kb with
    24 exons** — different genomic *scale*, so not a colinear copy.
  - Read-level: copy 3 = **0.1 reads/kb (14 reads/250 kb)** vs copy 1 = **10.1 reads/kb
    (682 reads)** — ~100× less coverage; copy 3's reads span up to 189 kb (read-through),
    not a dense gene fragmented by assembly.
  Borrowing copy 1's rich 24-exon/70 kb chain onto this sparse 250 kb locus would fabricate.

**Why borrow has no beneficiary (structural).** TOPO_BORROW assumes copies are homologous
with *identical structure* differing only in coverage. Real paralog families violate this:
structural/expression divergence is exactly what separates paralog copies. Copies similar
enough to borrow safely both assemble well (don't need it); copies thin enough to need it
are genuinely divergent or unexpressed (borrowing fabricates). This is *why* the `has_reads`
safety correctly yields 0 borrows on real data — it is not a missed opportunity. (Caveat:
two datasets, not a whole-genome sweep — but the structural argument generalizes, and these
cover the canonical multi-copy families: DAZ/RBMY/TSPY + GOLGA8/NBPF/TBC1D3/etc.)

**Verdict: retire TOPO_BORROW as an inert demonstrator.** The audit's #1 lever is resolved
(projection now correct) but has no real-data beneficiary. The genuine VG copy-resolution
value is elsewhere (copy-number-from-depth, confidence surfacing).

## Lever #3 BUILT: copy-number-from-depth — a diagnostic, not a hidden-copy oracle (2026-06-03)

Implemented `estimate_copies_from_depth` (pure, 6 unit tests) + `classify_family` wiring +
`family_depth_copies` GTF attribute, opt-in `RUSTLE_VG_DEPTH_COPYNUM` (default byte-identical,
commit d7d9123). Per-copy depth is attributed **multimap-aware** — a high-identity family's
bright copy lives mostly in `family.multimap_reads`, not `bundle.reads`, so each multimap read
is charged to its best-scoring copy (`reads_pc`) plus each copy's unique reads. Median
single-copy baseline (truncation-robust); optional `RUSTLE_VG_SINGLE_COPY_DEPTH` external unit.

**Honest validation (the build works; the SIGNAL is limited):**
- Equal-expression synthetic (1:1:1): **3.13** ✓.
- Real GOLGA8: largely **tracks the structural count** (FAM_12 2.00, FAM_14 3.89, FAM_6 4.87
  vs structural 2/4/4) and **flags gross imbalance** (FAM_16 55 vs structural 7).
- Collapse-detection (a true 2:1:1 where one bundle hides 2 copies): reads **3.20, not 4**.
  At 97% identity a copy's extra reads can't all be resolved to it (truncated 5′ reads miss
  diagnostic SNPs), so the estimate flattens toward the structural count.

**Why it can't reveal hidden copies (two fundamental limits):**
1. **Read-to-copy attribution is imperfect at realistic identity** — the very reads that
   would reveal a collapsed/dominant copy are the hardest to assign (this is the
   identifiability boundary the thesis is about). So depth-attribution converges toward the
   structural count instead of exceeding it.
2. **Depth conflates copy number with expression** — only equal-per-copy expression makes
   depth == copy number; real paralogs span ~100× expression.

**Verdict:** ship as an opt-in DIAGNOSTIC. Its real use is the **disagreement signal**:
`family_depth_copies` ≫ `family_n_copies` flags a dominant copy / possible collapse worth
human review (e.g. FAM_16). It is NOT a standalone hidden-copy counter. This is the same
pattern as TOPO_BORROW: clean in theory, bounded by the identifiability limit in practice —
an honest negative-leaning result, not a fabricated win.
