# Fuzzy Junction Skeleton Merging — Design

**Status**: DESIGNED, not yet implemented. Written 2026-09-15, following the chr20 ordinary-chromosome
assembler comparison (`bench/CHR20_ASSEMBLER_COMPARISON.md`), where our tool trailed StringTie on
gffcompare precision and most sensitivity levels. One suspected cause: `pass1_skeletons_robust_with`
(`src/rustle/vg_family/denovo_assemble.rs:256`) groups reads into "the same transcript" only on a
byte-exact intron-chain match (`BTreeMap<(&str, Vec<(u64,u64)>), ...>` keyed on exact `(start,end)`
tuples) — a few bp of real alignment jitter at a splice site is enough to split what should be one
transcript into several, each individually harder to match against the reference during scoring.

## Goal

Approach B from this session's brainstorming: **post-hoc skeleton merging**. After
`pass1_skeletons_robust_with` builds its exact-match skeletons exactly as it does today (unchanged), a
new, separate, opt-in step merges skeletons whose intron chains are structurally identical and differ
only by a small, pre-registered per-junction coordinate offset — before the merged result reaches the
`--gtf` writer.

**Why this shape, not a change to the grouping key itself (Approach A):** `pass1_skeletons_robust` is
called from 5 sites in `denovo_pipeline.rs` (lines 235, 2208, 2852, 3056, 4111), shared with the
multi-copy family / O1 detection pipeline this thesis's real, already-measured results are built on.
`family_detect.rs:540-546` documents this exact-match behavior as the KNOWN cause of real fragmentation
at real segmental-duplication loci (NOTCH2: 490 reads → 240 chains; SRGAP2C: biggest chain only 46 of
many more reads), already worked around downstream via `pick_locus_rep`'s `spliced_rep` — i.e. this exact
codebase already has precedent for reconciling fragmentation AFTER skeleton-building rather than
preventing it at the grouping key. Approach B follows that existing precedent and touches zero lines of
the shared, validated function; Approach A would not.

## Non-goal

Wiring this into the multi-copy family / O1 detection pipeline. This design targets only the general/de
novo `--gtf` emission path (no `--families` supplied). Extending it to the family pipeline, if ever
warranted, is separate future work requiring its own validation against the existing multi-copy results —
not decided here.

## Why an earlier fuzzy-matching attempt doesn't settle this

`docs/o1_ledger.md` §5f tested ±5bp junction snapping on one real 14-read NPIP-like locus and found no
difference (chains@exact == chains@±5bp), concluding the 14 distinct chains reflected genuine structural
conflict, not noise. A later finding (same document, ~line 3079) measured REAL alignment jitter on
simulated single-template reads and found it is much larger than 5bp (median donor-site gap 64bp; only
18.7% of gaps ≤5bp) — meaning §5f's ±5bp window was too narrow to fairly test the fuzzy-matching
hypothesis at all. That later finding explicitly left "junction placement tolerance at a realistic scale
(tens to low hundreds of bp, not 5)" untested on real reads. This design is that untested experiment, run
properly.

A separate, unrelated test (§6at) found ±20bp coordinate canonicalization net-worse for the
readthrough/mischain guard specifically (a different downstream filter, not this grouping key) — noted
here only as a reminder that "widen the tolerance" is not automatically safe, which is exactly why the
tolerance below is measured, not guessed.

## Architecture

```rust
/// Merge skeletons whose intron chains are structurally identical (same intron count, same strand) and
/// differ only by <= `tolerance_bp` at every corresponding donor/acceptor position. Single-linkage: if A
/// merges with B and B merges with C, all three end up in one group even if A and C alone exceed the
/// tolerance (same chaining pattern as this session's own `cluster_tie_partners`,
/// `src/rustle/vg_family/copy_discovery.rs`). Pure function over `Skeleton` values -- no BAM/env access,
/// unit-testable with synthetic skeletons.
pub fn merge_fuzzy_skeletons(skeletons: Vec<Skeleton>, tolerance_bp: u64) -> Vec<Skeleton>
```

Location: `src/rustle/vg_family/denovo_assemble.rs`, alongside `pass1_skeletons_robust` (needs direct
access to the `Skeleton` type defined there: `{chrom: String, start: u64, end: u64, n_reads: u32,
introns: Vec<(u64,u64)>, tied_seeded: bool, read_strand: Option<char>, ...}`, `denovo_assemble.rs:54-65`).

**Merge test** (both skeletons `a`, `b`):
1. `a.chrom == b.chrom`.
2. `a.introns.len() == b.introns.len()` (same number of introns -- never merge structurally different
   isoforms).
3. For every index `i`: `(a.introns[i].0 as i64 - b.introns[i].0 as i64).abs() <= tolerance_bp as i64` AND
   the same for `.1` (both donor and acceptor within tolerance, index-wise, since `introns` is already
   position-sorted).
4. Strand: `Skeleton.read_strand` is `Option<char>`, documented as "consulted ONLY for an unspliced model
   ... a spliced model's strand comes from its junction motifs" (`denovo_assemble.rs:61-64`) -- for a
   SPLICED skeleton (the only kind with a non-empty `introns` list, which is all this function ever
   considers) `read_strand` may be `None`. **Resolve during implementation**: confirm whether strand is
   available at the point this function is called (before or after `build_spliced_seq`'s motif-based
   strand assignment) and use whichever real strand signal is actually in scope there; if genuinely
   unavailable at this stage, drop the strand check from the merge test (a chain-count + per-junction
   tolerance match on the same chromosome, both directions, is unlikely to falsely conflate opposite-
   strand transcripts in practice, but do not silently skip this without one added unit test proving it
   using this codebase's own convention for such a caveat, matching `build_footprint_seq`'s documented
   placeholder note).

**Merge action**: union the two skeletons' reads (`n_reads = a.n_reads + b.n_reads`); for each intron
slot, keep the single most-common EXACT value among all pooled reads' own junction calls at that slot
(not an average or midpoint -- an invented coordinate that no read actually observed), ties broken by the
lowest coordinate (mirrors this codebase's existing tie-break convention, e.g. `build_footprint_seq`'s
`+`-on-tie strand rule and the discovered-copy strand majority vote's own tie rule). `start`/`end` and any
other `Skeleton` field: carry over from whichever input skeleton has more reads (ties: lower `start`) --
do not invent a new boundary-selection rule when this codebase already has an established quantile/robust
rule (`min_terminal_support`) for exactly this decision elsewhere; reuse that logic rather than average or
extremes if the implementation finds a clean way to do so, otherwise this simple more-reads/lower-tie rule
is an acceptable minimal choice for the merged skeleton's own boundary (its introns' identity is the part
this feature actually targets, not TSS/TES, which chr20's own already-tested `RUSTLE_TSS_SNAP`/
`RUSTLE_TES_EXTEND` showed no effect on -- see `bench/CHR20_ASSEMBLER_COMPARISON.md`'s follow-up section).

## Wiring (opt-in, general path only)

New env var `RUSTLE_JUNCTION_FUZZ_BP` (unset or `0` = off, matching this codebase's universal opt-in
convention -- e.g. `RUSTLE_TSS_SNAP`, `RUSTLE_TES_EXTEND`, `RUSTLE_JUNCTION_MAJORITY`). When set to a
positive integer, call `merge_fuzzy_skeletons(skeletons, value)` immediately after the ONE
`pass1_skeletons_robust` call site that actually serves `--gtf`'s pure de novo path (no `--families`
supplied) -- **resolve exactly which of the 5 call sites this is during implementation** by tracing which
one `detect_and_assign` reaches when `supplied` is `None` (per this session's own investigation,
`detect_and_assign` calls the whole skeleton-building front end unconditionally per
`denovo_pipeline.rs:2205-2248`, and short-circuits it entirely when a catalog IS supplied -- confirm the
exact line, do not guess). Do NOT wire this into the other 4 call sites without separately validating
against the multi-copy family test suite first -- out of scope for this design.

## Tolerance selection: pre-registration required, not a guessed constant

This is the part that directly answers the advisor's "did you tune this to pass" concern. Before any
merge code is evaluated against chr20's gffcompare/SQANTI3 numbers:

1. **Measure real jitter on chr20 itself** (not simulated data, not the earlier §5f/§3079 numbers, which
   were a different locus / simulated reads respectively): for every real chr20 read whose junction sits
   within some generous capture window (e.g. 500bp) of a real annotated RefSeq intron boundary, record the
   signed offset between the read's own CIGAR-derived junction coordinate and the annotated coordinate.
   Pool these into a real empirical distribution (a real, freshly-measured number, not assumed).
2. **Pick `tolerance_bp` by one fixed, stated rule, decided BEFORE step 3 runs**: the 90th percentile of
   the absolute offset distribution from step 1. (This specific percentile is chosen now, in this design
   doc, before any chr20 gffcompare number from this feature has been looked at -- do not change it after
   seeing whether it helps.)
3. Write the measured distribution, the resulting `tolerance_bp` value, and the exact command used to
   derive it to a new `docs/PREREG_junction_fuzz_2026-09-15.md`, committed to git BEFORE step 4 runs.
4. Only then: build `merge_fuzzy_skeletons`, wire it in behind `RUSTLE_JUNCTION_FUZZ_BP`, rerun the chr20
   `bakeoff_chr20_ours.sh` with the flag set to the pre-registered value, rescore with gffcompare/SQANTI3
   exactly as `bench/CHR20_ASSEMBLER_COMPARISON.md` already does, and report the result HONESTLY --
   including if it does not help, following this project's own standing culture (§6l6's own acceptance
   test reported a negative result plainly; this feature does the same if that is what happens).

## Testing plan

1. **Unit tests for `merge_fuzzy_skeletons`** (synthetic `Skeleton` values, no BAM/reads needed):
   - Two skeletons with the same intron count, every junction within `tolerance_bp` -> merge into one,
     `n_reads` summed.
   - Two skeletons with a DIFFERENT intron count -> never merge, regardless of how close the shared
     junctions are.
   - Two skeletons with the same intron count but ANY one junction beyond `tolerance_bp` -> never merge.
   - Three skeletons A-B-C where A-B and B-C are each within tolerance but A-C alone is not -> all three
     merge into one group (single-linkage chaining, matching `cluster_tie_partners`'s own tested
     behavior).
   - `tolerance_bp = 0` -> behaves identically to no merging at all (regression guard: this feature must
     never change output when explicitly set to zero, mirroring every other opt-in flag's off-state
     contract in this codebase).
2. **Byte-identical-when-unset**: `RUSTLE_JUNCTION_FUZZ_BP` unset must produce byte-identical `--gtf`
   output to the current binary -- add a regression test in the same style as
   `discover_copies_off_by_default_is_byte_identical` (`tests/copy_assign_families.rs`).
3. **Real-data acceptance**: after the PREREG value is measured and committed, rerun the chr20 comparison
   with the flag on, report the real gffcompare/SQANTI3 delta in `bench/CHR20_ASSEMBLER_COMPARISON.md` as
   a new dated section, honestly, whichever way it goes.

## Global constraints

- `RUSTLE_JUNCTION_FUZZ_BP` unset (or `0`) must produce byte-identical output to the current binary.
- Never wired into the multi-copy family / O1 detection pipeline in this design -- general/`--gtf` de
  novo path only.
- The tolerance value is measured from real chr20 data and fixed BEFORE its effect on the metric is
  observed -- no re-deriving a different tolerance after seeing a disappointing first result without
  saying so explicitly and re-registering.
- No genome-only discovery / no arbitrary new threshold without a stated reason -- the tolerance's
  provenance (90th percentile of measured real jitter) is that stated reason.
- Reuse existing tie-break/majority conventions (most-common exact value, ties to lowest coordinate) for
  the merged junction position -- do not invent a new averaging or interpolation rule.
