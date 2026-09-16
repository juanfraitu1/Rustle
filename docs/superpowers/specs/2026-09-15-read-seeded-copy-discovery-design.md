# Read-Seeded Copy Discovery — Design

**Status**: DESIGNED, not yet implemented. Written 2026-09-15, following §6l5's finding
(`docs/o1_ledger.md`) that 90% of real AS-tied reads in the audited catalog tie against a real BAM
placement that sits outside every catalog copy, not against another catalog member.

## Goal

Phase 1 of "build the assembler side so it shines on multi-copy families and incomplete paralogs"
(user request, this session). Give `copy_assign --families` a way to surface candidate copies that its
own AS-tie evidence already points at but the supplied catalog is missing — closing exactly the gap
§6l5 measured, not a general de novo assembly rewrite.

**Non-goal (deferred, Phase 2):** porting `bench/dna_sd_atoms.py`'s DNA self-alignment ("SD atom")
node-discovery into production. That is a separate, larger design or write-up on its own once Phase 1
ships and its acceptance test (below) is measured on real data.

**Non-goal:** matching StringTie/FLAIR/IsoSeq's general isoform-recall bar (sub-project A from this
session's brainstorming). Out of scope for this design entirely.

## Why the two existing "rescue" mechanisms don't already solve this

Investigated this session (fork report, not re-quoted in full — see chat history / re-derive if needed):
`rescue_thin_loci_iterative` (`src/rustle/vg_family/rescue_pipeline.rs`) and `--tied-seed`'s `TSFAM`
existence-only append (`src/rustle/vg_family/denovo_pipeline.rs:1230-1237`) are both real, tested code —
but both live inside `detect_and_assign` (`denovo_pipeline.rs`), which is `copy_assign`'s **own internal
de novo detection path**, used only when no `--families` catalog is supplied. Every real catalog audited
this session was built externally by `gw_family_catalog` and fed in via `--families`, which **skips
`detect_and_assign` entirely** (`denovo_pipeline.rs:2429-2460`) — so for the actual, currently-used
pipeline, both mechanisms are unreachable, not merely under-tuned. Even if reachable, `thin_loci`
(`rescue_pipeline.rs:57-120`) seeds from a read's own **spliced, multi-exon** primary interval — it has
no concept of an AS-tie partner position, so it would not have found §6l5's candidates even if wired in.

This design therefore adds a **new** mechanism inside `copy_assign --families`'s own code path, where
AS-tie evidence is already computed, rather than trying to repair either existing one.

## Architecture

One new opt-in CLI flag, `--discover-copies` (bool, default off — byte-identical output when unset, same
convention as every other opt-in flag in this binary). When set, after the normal per-region assignment
pass:

1. For every read already flagged AS-tied (`as_margin == 0.0` in the existing `AsEvidence`, computed by
   `as_evidence_per_read`, `src/bin/copy_assign.rs:1613`), re-examine that read's **full placement set**
   (every primary/secondary `BamRead` sharing its name within the region, the same population
   `as_evidence_per_read` already groups via its internal `by_name` map, `copy_assign.rs` — see that
   function's body for the exact grouping key/filter, `exclude_supplementary`).
2. For each placement scoring at the read's max AS (there may be more than 2 — do not assume exactly
   two), test whether its `(chrom, ref_start, ref_end)` falls inside any `CatalogCopy` span
   (`src/rustle/vg_family/catalog_input.rs:51-73`) belonging to the read's own family (from the already-
   loaded `RegionFamilies`/`CatalogFamily.copies`, `copy_assign.rs:1278`, `1419` `load_supplied_families`).
   A placement with NO containing copy is an **out-of-catalog tie partner**.
3. Cluster out-of-catalog tie-partner positions across all reads in a family: two positions merge into
   one candidate site if their genomic intervals overlap or sit within 500 bp of each other (a fixed,
   documented constant — not tuned on this data; revisit if the acceptance test below shows clusters
   splitting a single real locus). Each resulting cluster's supporting-read count = number of **distinct
   read names** contributing a placement to it.
4. Keep clusters with `n_supporting_reads >= 2` — the same bar `read_supported_columns` already uses
   for PSV columns ("keeps a column only if reads show ≥2 alleles at ≥2 reads each", §6aj) — as the
   admission threshold. This is a reused project convention, not a new tuned constant.
5. Emit `<out>.discovered_copies.tsv`: one row per kept cluster —
   `family_id, chrom, cand_start, cand_end, n_supporting_reads, read_names (comma-joined),
   nearest_catalog_copy_tid, nearest_catalog_copy_distance_bp`. This file is a REPORT, not a live catalog
   edit — `copy_assign`'s own assignment output for this run is unaffected by anything discovered (matches
   the existing `--tied-seed` convention: discovery/rescue never perturbs the primary run's byte-for-byte
   output).

## Two-pass workflow (explicit, not automatic)

`copy_assign --families cat.tsv ... --discover-copies` → inspect `<out>.discovered_copies.tsv` → for
accepted clusters, append a new row to `cat.tsv` (family_id, next copy_idx, a synthetic tid, chrom,
cand_start, cand_end, strand `?` pending a real strand call — see Open Question below, n_reads =
n_supporting_reads, exons = single-block `[cand_start, cand_end)`, core_hull/locus = `NA`, partner =
`false`) → re-run `copy_assign --families cat.tsv` (no `--discover-copies`) normally with the enlarged
catalog. No code changes needed for the merge step in Phase 1 — a human or a small follow-up script
(`bench/merge_discovered_copies.py`, not specified further here) does the append; this design only
specifies the discovery half.

## Data structures (new, `src/rustle/vg_family/copy_discovery.rs`)

```rust
/// One candidate copy inferred from AS-tied reads whose best-scoring placement has no catalog home.
#[derive(Clone, Debug)]
pub struct DiscoveredCopy {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_supporting_reads: usize,
    pub read_names: Vec<String>,
    pub nearest_copy_tid: String,
    pub nearest_copy_distance: u64,
}

/// Cluster out-of-catalog AS-tie partners into candidate copies for one family.
/// `tied_reads`: (read_name, Vec<(chrom, ref_start, ref_end)>) — every max-AS placement of every
/// AS-tied read in this family, already filtered to exclude any placement that DOES fall inside a
/// `family.copies` span (that filtering happens in the caller, not this function — keep this function a
/// pure clustering primitive, unit-testable without BAM/catalog plumbing).
pub fn cluster_tie_partners(
    tied_reads: &[(String, Vec<(String, u64, u64)>)],
    family: &CatalogFamily,
    merge_distance: u64,
    min_support: usize,
) -> Vec<DiscoveredCopy>
```

`copy_assign.rs` gets one new function, `discover_copies_for_family`, that (a) filters the region's tied
reads down to out-of-catalog placements per family and (b) calls `cluster_tie_partners`; and one new
writer, mirroring the existing `<out>.X.tsv` writers already in `main()`, gated on `args.discover_copies`.

## Open question to resolve during implementation (not blocking the plan, flag and pick one)

Strand for a discovered copy: a raw AS-tied placement's own SAM strand flag is available per-record, but
may disagree across the supporting reads (some `+`, some `-`, e.g. antisense contamination). Simplest
correct rule: majority vote across supporting reads' own `reverse` flags (matches `majority_read_strand`'s
existing convention elsewhere in this codebase per `denovo_pipeline.rs:3957-3961`'s own comment); ties
default to `+` with the same historical-placeholder caveat already documented for footprints
(`build_footprint_seq`, `denovo_assemble.rs`). Use that rule; do not invent a new one.

## Acceptance test

Re-run on the real §6l5 substrate: `arm_f2`'s 6 low-identity families (GWFAM55, GWFAM66, GWFAM96,
GWFAM104, GWFAM113, GWFAM118), `npip_cat/npip3.bam`, catalog `mec/psv.*`'s own `copies.tsv` (or the
equivalent `arm_f2/cat.copies.tsv` — confirm which one `mec/psv.*` was actually built from before running;
this session's own investigations used both names for what may or may not be the same file, re-derive
rather than assume). Expect `<out>.discovered_copies.tsv` to surface a cluster at or near
`NC_073242.2:21,674,468` (GWFAM55, the §6l5-verified example, 31.7 kb before copy 2) with
`n_supporting_reads >= 2`. If it does not appear, the merge-distance or min-support constants need
revisiting, not a redesign.

Do NOT quote a "coverage improved" headline from this alone — the acceptance test is existence-only
(does the known real candidate get discovered), not a recall/precision claim over the whole genome. A
broader validation (how many of the 65 out-of-catalog tie-partner reads from §6l5 get explained after one
merge-and-rerun cycle) is a natural follow-up once Phase 1 ships, not part of this design's acceptance bar.

## Testing plan

1. Unit tests for `cluster_tie_partners` (`src/rustle/vg_family/copy_discovery.rs`'s own `#[cfg(test)]`
   module): synthetic tied-read position lists — merges positions within `merge_distance`, keeps distinct
   clusters beyond it, respects `min_support`, correctly excludes positions already inside a synthetic
   `CatalogFamily`'s copy spans (i.e., confirm the CALLER's filtering contract holds if a caller mistake
   passes an in-catalog position — decide whether the function should also defensively re-check this
   itself; recommend yes, cheap and removes a whole class of caller bugs).
2. Integration: run `copy_assign --families --discover-copies` on the real GWFAM55 substrate above,
   assert the known cluster appears with the right support count (independently re-derivable via the same
   `samtools view <bam> <region>` method already used to verify it by hand this session).
3. Regression: run the full existing `--families` test suite (`tests/copy_assign_families.rs` and
   friends) with `--discover-copies` OFF (the default) and confirm byte-identical output to before this
   change — this flag must never alter the primary assignment path.

## Global constraints

- Default OFF; `--discover-copies` unset must produce byte-identical output to the current binary.
- Never mutates the input catalog or the current run's own assignments — report file only (two-pass,
  per user's explicit choice this session over an automatic single-pass merge).
- Reuse the `>=2` independent-read-support convention already established for PSV columns; do not invent
  a new threshold without a stated reason.
- No genome-only discovery: every candidate must be seeded by a real AS-tied read's own real alignment
  record — never a blind genome scan (matches the project's standing `NO genome-only discovery` scope
  rule, `docs/PREREG_core_definition_2026-09-12.md` era decisions).
