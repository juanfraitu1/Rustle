# O3 Flag-Pass Integration Design

**Status**: design, not yet implemented. Written 2026-09-10.

## Goal

O3 (THESIS_OBJECTIVES.md: "Detect + flag expressed transcript paths not explained by represented
reference copies, STRATIFIED by whether the orphaned reads have anywhere to go") currently exists only
as a real, measured, but disconnected mechanism: three standalone Python scripts
(`bench/o3_flag_pass.py`, `bench/o3_cut_certificate.py`, `bench/o3_reconstruct.py`) glued to
`copy_assign`'s own output by a directory convention (one `fam_*` directory per family, each holding a
prior `copy_assign` run's `copies.tsv` + `A.assignments.tsv`). There is no single command; running O3
today means: run `copy_assign` once per family into a sweep directory, then run `o3_flag_pass.py` over
that directory.

This design ports the **flag-pass detector only** (`bench/o3_flag_pass.py`; the cut-certificate and
reconstruction scripts are explicitly out of scope, per user decision) natively into the Rust pipeline,
as a new opt-in flag on `copy_assign`, so that:
1. O3 runs as part of the same command that already produces O1/O2's output (**one flow**).
2. The detector logic lives in Rust, not Python, reusing data `copy_assign` already holds in memory
   instead of round-tripping through disk (**ported to Rust**, and avoids the "two independent
   implementations can drift" risk already flagged elsewhere in this codebase, e.g. E3's "the lift rule
   exists twice").
3. Its output lands on **existing** `copy_assign` output files wherever a natural home exists, rather
   than only ever appearing in a disconnected sidecar (**surfaced in existing outputs**).

### Caveat: expected flag rate on the matched substrate is near-zero, and that is correct, not a bug

The standard gorilla substrate (fibroblast IsoSeq) is proven same-CELL-LINE as the mGorGor1 assembly
itself (SRA-accession match + a 13.5× hom-alt-vs-het internal control, `project_o3_matched_individual.md`).
On this substrate, a reference-absent copy has very little room to exist by construction — the reference
*is* this individual's genome. Measured, not assumed: the shipped RNA-side detector already fires
**0/915** on it, and an independent DNA-side assembly-vs-assembly check (the field's own S1 standard)
found **0/817** collapse-shaped deficits, a result the literature's own prior (Yoo/Rhie 2025: "1–2 Mbp
collapse per haplotype") predicts in advance (0.47–0.94 expected in that compartment). **Do not read a
near-zero `missing_copy` flag rate from `--flag-missing-copies` on this substrate as evidence the port
is broken — it is the expected result of running a reference-absence detector on a reference that
already contains the tested individual.**

The one place this project found a real, if thin, residual even under exact-match conditions is a
one-time exploratory pass, not this detector: 3 candidate loci at ~0.94 identity to every haplotype of
the diploid assembly, not explained by repeat content, not rescued by either parental haplotype — but
with only n=3, no PCR validation, and a live undisambiguated confound (cultured fibroblast lines
accumulate somatic drift from the sequenced germline DNA independent of any real "missing copy").
TSPY (Rhie et al. 2023, *Nature* 621) is the concrete proof that the large-magnitude version of this
phenomenon — 10–40 copies across a human population against one reference individual — needs
**non-matched-individual** data to appear at all. See the "Detecting the real signal" addendum below,
which is a second, complementary tool this design does not itself build.

## Non-goals

- Porting `o3_cut_certificate.py` (scores O1 MCL cuts) or `o3_reconstruct.py` (patches consensus at
  flagged sites) — both stay separate Python tools, run by hand when needed.
- Resolving D5 (RNA admission of unannotated loci into O1) or D6's four sub-decisions
  (`docs/OPEN_ITEMS_2026-09-09.md`). This design surfaces the same orphan-read stratification the
  Python script already produces (which feeds those decisions) but does not decide them. Nothing here
  changes O1's catalog or O2's assignment.
- Changing the existing `--absent-copies` flag (`src/rustle/vg_family/absent_copy.rs`, "gate 5") — a
  different, older, unrelated reference-absent reconstruction attempt, already default-off and never
  used in production per `docs/O3_STATUS.md`. This design does not touch it, and the new flag uses a
  clearly distinct name so the two are never confused.
- Changing anything about the default (flag-unset) output of `copy_assign` in any way. Every new field
  in this design is additive and gated.

## Architecture

### Why not put this directly in `copy_assign.rs`

`copy_assign.rs` is already 4,600+ lines. The detector's logic (build a local sequence window for a
candidate locus, realign a batch of reads to it via minimap2, count consistent mismatch sites, run a
Poisson tail test) is a self-contained unit with a clean interface, exactly the shape of the existing
`src/rustle/vg_family/absent_copy.rs` and `linearize.rs` modules that `copy_assign.rs` already calls
into rather than inlines. It goes in a new module for the same reason.

### Why not a separate binary

A separate binary reading `copy_assign`'s completed TSV output would still need its own `--bam`/`--fasta`
access and would still have to re-fetch and realign the same reads from the BAM — it does not save real
work, and it reintroduces exactly the "two commands glued by files" seam this design exists to remove.

### The two-phase split (why it's required, not a stylistic choice)

The Python detector's final flag is **not** a fixed per-pair threshold: `flags.tsv`'s aggregate step
computes `n_pairs` = the total number of testable pairs across the **entire run** (every family, every
contig), then `thr = alpha / n_pairs`, and only then labels each pair `missing_copy` (a Bonferroni
correction across everything tested this run). This means:

- **Phase 1 — per-family raw computation** can run inside `copy_assign`'s existing per-region parallel
  `compute()` closure (same closure that already builds each region's `fams`/`transcripts`/certificate).
  It needs no data from other families. Produces one `RawPair` per (family, candidate-copy-Y) with
  ≥3 origin-rejected reads, carrying the **uncorrected** p-value, plus `OrphanLocus` candidates.
- **Phase 2 — genome-wide aggregation** must wait until every region has been computed. `copy_assign`
  already has exactly this shape of dependency and exactly this solution: the `productive`/`orf_aa` GTF
  attribute is "RELATIVE to the family's best ORF, which is only known once every region has been
  drained... stamped here, in a second pass over the finished GTF lines" (existing comment, `copy_assign.rs`
  around the productivity pass). O3's aggregation is a second such pass, run at the same point (after
  the serial `for work in works.iter()` drain, before final output is written): collect every `RawPair`
  from every region, count `n_pairs`, compute the threshold, assign final flags.

### Data already available vs. new data needed

Already resident per region, reusable as-is:
- `verdict: HashMap<&str, &AssignRow>` — `AssignRow` (line ~1604) already carries `origin_rejected: bool`,
  `n_candidates: usize` (0 ⇒ orphan, matching the Python's `n_candidates=='0'` check), `catalog_copy_idx:
  String`, `status: &'static str` — everything the Python derives from `A.assignments.tsv` is already an
  in-memory struct field, keyed by read name exactly like the existing B2/read-provenance classify loop
  already does.
- `fams[fw].copy_spans` — each family's own copies' `(chrom, start, end)`, matching the Python's `cp`
  dict built from `copies.tsv`.
- `genome_for(contig)` — the read-only cached genome accessor `copy_assign.rs` already uses for exon
  sequence fetches, reusable for building Y's local window FASTA.

New, needed once per run (built before the sweep starts, from data already parsed):
- `all_units_by_chrom: BTreeMap<String, Vec<UnitRef>>` — every family's copies, flattened across the
  whole `region_families` map (already fully loaded by `load_supplied_families` before the sweep begins,
  at `copy_assign.rs:1799`), re-indexed by chrom. Needed for the `OtherFamily` classification of a
  candidate orphan locus (does it overlap a DIFFERENT family's unit) — this is inherently cross-family,
  so it cannot come from any single region's own `fams`.
- `genes_by_chrom: BTreeMap<String, Vec<(u64,u64)>>` — reuses the existing `parse_annotation()` (already
  parses `--gff`/`--bed` into `Vec<(String,u64,u64)>`) grouped by chrom, for the `AnnotatedNoUnit` vs
  `Unannotated` classification. **Scope note**: `parse_annotation` returns no gene name/biotype today;
  the initial port reports orphan-locus overlap as a count (`n_genes_overlapping`), not named genes —
  matching the Python's own display fields is a follow-up, not required for the flag itself.

Both are built once, shared read-only across every parallel region worker (same pattern as the existing
`genome_cache`/`bam_cache`).

### The `family_join.tsv` sequencing problem (concrete, not hand-waved)

`copy_assign.rs`'s `family_join.tsv` rows (`join_rows: Vec<String>`, pushed at line ~2590 during the
serial `for work in works.iter()` drain) are pre-formatted, already-joined strings, written directly at
the end (line ~3786). O3's final flag is only known from the Phase-2 aggregation, which itself only runs
after that same serial drain completes — so a `RawPair`'s flag cannot be known at the moment its
`family_join.tsv` row string is built.

Resolution: change `join_rows` from `Vec<String>` to a small struct carrying the pre-formatted line
alongside its join key (`family_id`, `catalog_copy_idx`) — e.g. `JoinRow { line: String, family_id:
String, copy_idx: String }` — so the write loop (line ~3786) can look up this row's finalized O3 flag (by
`(family_id, copy_idx)` in a `HashMap` built from Phase 2's output) and append the O3 columns to `line`
at write time, without changing anything about how the existing columns are built. This is a small,
mechanical, purely-additive refactor of `join_rows`'s element type — called out explicitly here because
it is exactly the kind of detail that's easy to hand-wave in a design and then discover mid-implementation.

## Components

### New module: `src/rustle/vg_family/o3_flag_pass.rs`

```
pub enum Class { Divergent, Structural }
pub enum LocusClass { OtherFamily, AnnotatedNoUnit, Unannotated }
pub enum Flag { MissingCopy, Untestable, None }

pub struct RawPair {
    pub family_id: String,
    pub copy_idx: String,          // matches AssignRow::catalog_copy_idx's namespace
    pub is_partner: bool,          // CatalogCopy::partner (the parsed form of the raw member_status
                                    // column; "NA"-vs-real-string granularity is not preserved by the
                                    // existing parser and is not needed for the flag decision itself)
    pub n_rejected: usize,
    pub n_aligned: usize,
    pub covered_kb: f64,
    pub n_sites: usize,
    pub ctl_n: usize,
    pub ctl_covered_kb: f64,
    pub ctl_n_sites: usize,
    pub p_uncorrected: Option<f64>,  // None when kb or ctl_kb is 0 (Python's float('nan') case)
    pub class: Class,
    pub med_mismatch: f64,   // median per-read mismatch count, over the REJECTED reads only
    pub med_unaligned: f64,  // median per-read unaligned-base count, over the REJECTED reads only
}
// Class rule (mirrors the Python exactly): Divergent when med_mismatch > med_unaligned (substitution
// pattern dominates), else Structural (unaligned/indel pattern dominates, or med_mismatch is NaN
// because no rejected read realigned at all).

pub struct OrphanLocus {
    pub chrom: String, pub start: u64, pub end: u64,
    pub n_reads: usize, pub n_orphans: usize,
    pub class: LocusClass,
    pub n_genes_overlapping: usize,
    pub other_family_units: Vec<String>,  // "family_id:copy_idx", capped at 3 (Python also embeds the raw
                                       // member_status string here; dropped, see is_partner note above)
}
// LocusClass precedence (mirrors the Python exactly): OtherFamily if the locus overlaps ANY other
// family's unit, else AnnotatedNoUnit if it overlaps a `--gff` gene/pseudogene, else Unannotated.
// Checked in that order — a locus overlapping both another family's unit AND a gene is OtherFamily,
// never AnnotatedNoUnit.

pub struct FlaggedPair { pub pair: RawPair, pub flag: Flag, pub p_corrected_threshold: f64 }

pub struct O3Params { pub alpha: f64, pub max_reads: usize, pub min_reads: usize } // defaults 0.001/500/3

/// Phase 1. Runs inside one region's compute(): groups this family's origin-rejected reads by best
/// candidate copy Y (≥ params.min_reads), realigns them + Y's own certificate-accepted reads (control)
/// to Y's local window via minimap2 -x splice (mirrors minimap2_msa_pair's temp-file/nonce pattern,
/// `copy_assign_pipeline.rs:131-169`), computes consistent-mismatch-sites/kb for both. `p_uncorrected`
/// is computed HERE (not in Phase 2 — Phase 2 only decides the acceptance threshold): rate =
/// `max(ctl_n_sites, 1) / ctl_covered_kb` (floored at 1 site so a control with zero observed sites
/// never claims a zero rate, which would trivially fail every test), `lambda = rate * covered_kb`,
/// `p_uncorrected = poisson_tail(n_sites, lambda)`, `None` when `covered_kb` or `ctl_covered_kb` is 0.
/// Also scans this family's rejected/orphan reads whose primary lands outside every one of the
/// family's own units (±200kb, matching the Python's window) into read clusters, classified via
/// `all_units_by_chrom`/`genes_by_chrom` (precedence: see `LocusClass` above).
pub fn compute_family_raw(
    family_id: &str,
    copy_spans: &[(String, u64, u64)],   // fams[fw].copy_spans, chrom/start/end per copy_idx
    is_partner: &[bool],                  // CatalogCopy::partner, parallel to copy_spans
    bam_reads: &[BamRead],
    verdict: &std::collections::HashMap<&str, &AssignRow>,
    genome: &GenomeIndex,
    all_units_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>>, // start,end,family_id,copy_idx
    genes_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
    params: &O3Params,
) -> anyhow::Result<(Vec<RawPair>, Vec<OrphanLocus>)>;

/// Phase 2. Runs once after every region has drained: counts n_pairs across every RawPair with
/// p_uncorrected.is_some() from the WHOLE run, threshold = alpha / n_pairs, labels each pair. Pure
/// function (no I/O), unit-testable directly against hand-computed reference cases.
pub fn finalize_flags(all_pairs: &[RawPair], alpha: f64) -> Vec<FlaggedPair>;
```

A new `poisson_tail(k: usize, lam: f64) -> f64` helper (P(X ≥ k) for X ~ Poisson(lam), same closed-form
sum the Python uses) lives in this module; it reuses the existing hand-rolled `lgamma` in
`allele_specific_junctions.rs` if that function is made `pub(crate)` (currently private to that module —
a one-line visibility change, or a duplicate ~10-line implementation if keeping modules fully
decoupled is preferred; either way, unit-tested against a handful of hand-computed values matching the
Python's own `poisson_tail`).

### `copy_assign.rs` changes

- New CLI flags: `--flag-missing-copies` (bool, default false), `--o3-alpha <f64>` (default 0.001),
  `--o3-max-reads <usize>` (default 500) — mirrors the Python's own `--alpha`/`--max-reads` tunables
  rather than hardcoding them, matching this codebase's convention of exposing measured constants.
- `all_units_by_chrom`/`genes_by_chrom` built once, right after `load_supplied_families` returns
  (line ~1799), only when `--flag-missing-copies` is set (zero cost otherwise).
- `RegionWork` gains two fields: `o3_raw_pairs: Vec<RawPair>`, `o3_orphan_loci: Vec<OrphanLocus>` — empty
  `Vec`s when the flag is unset (byte-identical struct size cost is negligible and matches how
  `transcripts: Vec<TranscriptRec>` is already "empty unless `--gtf`").
  The `compute()` closure calls `o3_flag_pass::compute_family_raw(...)` once per family in the region,
  only when `args.flag_missing_copies`.
- Serial drain (`for work in works.iter()`): accumulate `o3_raw_pairs`/`o3_orphan_loci` into two
  run-level `Vec`s, the same way `famcn_rows`/`family_rows` already accumulate.
- After the drain loop, before any output is written: if the flag is set, call
  `o3_flag_pass::finalize_flags(&all_raw_pairs, args.o3_alpha)`, build a
  `HashMap<(family_id, copy_idx), FlaggedPair>` from the result.
- `join_rows` becomes `Vec<JoinRow>` (see sequencing section above); the `family_join.tsv` write loop
  appends `\to3_flag\to3_class\to3_rate_per_kb\to3_p\to3_n_rejected` (tab-prefixed, only when the flag is
  set — the header line itself only gains these columns when set, so the unset case is the exact
  existing 9-column file).
- New `<out>.o3_candidate_loci.tsv`, written only when the flag is set: header
  `chrom\tstart\tend\tn_reads\tn_orphans\tclass\tn_genes_overlapping\tother_family_units`, one row per
  accumulated `OrphanLocus` — same information as the Python's `loci.tsv`, in a format any existing
  `bench/` consumer of that shape can read unchanged.

## Error handling

- **< min_reads (3) rejected reads for a candidate Y**: skip that pair entirely (not emitted as
  `Untestable` — matches the Python's `len(names) < 3: continue`).
- **Control set < min_reads**: still test (Python's fallback: `ctl_n = len(ctl_names)`, `ctl_kb = 0.0`,
  `ctl_sites = {}`), which forces `p_uncorrected = None` downstream (both `kb` and `ctl_kb` must be
  nonzero for a real p-value) — becomes `Flag::Untestable` in Phase 2, exactly matching the Python's
  `untestable` class.
- **`minimap2` subprocess failure for one Y**: log a warning (`[o3-flag-pass] realignment failed for
  {family_id}:{copy_idx}: {err}`) and skip that pair (contributes neither a flag nor an error to the
  whole run) — O3 is best-effort detection layered on top of an already-complete O1/O2 result; one
  failed realignment must never abort a multi-hour genome-wide sweep. This is a deliberate asymmetry
  from `copy_assign`'s treatment of, e.g., a malformed `--families` file (which does abort, at the
  boundary, before any work starts) — O3's realignment failures are a per-pair runtime condition, not an
  input-validation failure.
- **Temp file collisions under `--region-threads > 1`**: reuses the exact pid+atomic-nonce naming
  already established in `minimap2_msa_pair` (`copy_assign_pipeline.rs:139-143`) — no new risk class.
- **`--flag-missing-copies` without `--families`**: the detector has no catalog copies to test against;
  bail with a clear error at argument-validation time (same point `--gtf-copy-set` and other
  `--families`-dependent flags are already checked), rather than silently emitting nothing.

## Testing

- **Unit tests** (in `o3_flag_pass.rs`, following this file's existing test-module convention):
  `poisson_tail` against hand-computed values (including edge cases `k=0`, `lam=0`); `finalize_flags`
  against a small synthetic `Vec<RawPair>` where the correct threshold and resulting flags can be
  computed by hand; the `OrphanLocus` classifier (`OtherFamily` > `AnnotatedNoUnit` > `Unannotated`
  precedence, matching the Python's `if other else (... if g else ...)`) against constructed inputs for
  each of the three classes.
- **Byte-identity gate**: `--flag-missing-copies` unset ⇒ `copy_assign`'s entire output (GTF,
  `assignments.tsv`, `family_join.tsv`, everything) byte-identical to today, verified the same way every
  other flag this session was (build → run with/without the flag on the human chr16 substrate → diff).
- **Reproduction gate against the Python's own numbers**: re-run `bench/o3_flag_pass.py` once more on
  whichever substrate the ledger's §6fm-§6ft numbers came from (or the closest available equivalent),
  then run the new `--flag-missing-copies` on the identical input, and diff: pair count, `class`
  (divergent/structural) per pair, and final `flag` label should match; `p_uncorrected` should match to a
  stated float tolerance (independent implementations of the same closed-form sum). Any mismatch is
  investigated as a genuine discrepancy before this ships, not waved through.
- **Real-substrate smoke test**: run on the gorilla MCL1 substrate and the held-out human chr16 hard
  locus already established this session, confirm the new `family_join.tsv` columns and
  `o3_candidate_loci.tsv` are well-formed and the run completes in reasonable time (the realignment step
  adds real minimap2 subprocess cost per candidate Y — worth a rough wall-clock number in the plan, not
  assumed free).

## Follow-on: cross-individual differential (detecting the TSPY-style signal)

Out of scope for the flag-pass port itself, but reuses it entirely — recorded here rather than in a
separate spec because it is one more invocation of the same module, not a new detector.

### Why this is not speculative — the data already exists

Two independent, already-built gorilla substrates exist on disk today:
- `/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam` — **matched**: proven (SRA
  accession match + a 13.5× hom-alt-vs-het internal control) to be the exact cell line the mGorGor1
  assembly was built from.
- `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_ds.bam` (parent: `GGO_mm.bam`) — **OR6737, testis, a
  different animal**. This is not a side dataset: it is the substrate behind `gw_units_v3`, the
  genome-wide catalog this project's O1 work (including this session's own C6 hold-out pick) is built
  from. A depth-matched fibroblast replicate, `o1_replicate/fibro_ds.bam`, already exists for exactly
  this kind of paired comparison.

A prior attempt at exactly this comparison (ledger §4l, "cross-substrate replication," 2026-08-23) found
O1's family relation transports well across the pair (87.06% edge recovery, ARI 0.9707), but its
per-family **read-presence** screen (does a copy have ≥1 read in each tissue) came back null with no
directional bias (33/52 differ, 16/17 split, binomial p=1.0) — and the ledger's own verdict is explicit:
*"the read-presence screen cannot test the hypothesis — it measures expression."* A copy genuinely absent
from one individual and a copy merely silenced in one tissue look identical under presence/absence. No
follow-up built a copy-number-*shaped* test on this pair — that gap is what this section closes.

### The design: diff two flag-pass runs, not two read counts

Run `copy_assign --families mcl_ann/gw_units_v3.units.tsv --flag-missing-copies` twice against the SAME
catalog — once with `GGO_ds.bam` (testis/OR6737), once with `fibro_ds.bam` (matched fibroblast, already
depth-matched) — then compare each copy `Y`'s flag between the two runs instead of comparing raw read
counts. This is the fix for exactly what killed the presence screen: `Flag::MissingCopy` already tests
the *shape* of a mismatch pattern (consistent-site density against Y's own control), not raw depth, so a
tissue-driven expression difference that affects both the test and control populations similarly should
not by itself flip the flag, while a genuine between-individual copy-count difference plausibly would.

Comparison categories (deliberately distinct from a simple boolean diff, to avoid re-making the §4l
mistake of conflating "no data" with "no signal"):

| Testis flag | Fibroblast flag | Category | Interpretation |
|---|---|---|---|
| `MissingCopy` | `None` | **Candidate differential** | Both arms had adequate testable data; only one shows the divergence signature — the interesting case. |
| `None` | `MissingCopy` | **Candidate differential** | Same, other direction. |
| `MissingCopy` | `MissingCopy` | Shared | Signal in both individuals against this reference — an O1 catalog-completeness question, not a between-individual one. |
| `MissingCopy` or `None` | `Untestable` (or vice versa) | Inconclusive | One arm lacked enough control/rejected reads to test at all; report separately, never counted as a confirmed differential. |
| `Untestable` | `Untestable` | No information | Drop. |
| `None` | `None` | No signal | Not flagged. |

Only the **Candidate differential** row is the TSPY-style hypothesis; everything else is reported for
transparency but not claimed as a finding.

### Scope: a comparison script, not new detector code

The diff itself is a join over two already-written `<out>.family_join.tsv` files on `(family_id,
copy_idx)`, comparing the `o3_flag` column — no realignment, no statistics of its own. This is exactly
the shape of the existing `bench/` comparison scripts (e.g. `isoform_bakeoff.py` diffing multiple tools'
GTFs), so it belongs there as a small Python script (`bench/o3_cross_individual_diff.py`), not as more
Rust: the detector logic already lives in Rust (this design), the cross-run comparison is a lightweight
join that gains nothing from a native port.

### What's NOT established, and must not be overclaimed

- **This is exploratory, not validated.** There is no known-true-positive gorilla case (no gorilla TSPY
  equivalent with an independently confirmed answer) to check the comparison against — unlike the
  human paper's TSPY, where ddPCR and AmpliCoNE independently confirm the count. Any "candidate
  differential" result is a hypothesis for follow-up (targeted PCR, a third individual, or literature
  cross-reference for that specific gene family), not a proven missing copy.
- **The catalog-transport assumption needs its own sanity check before trusting any result**: `gw_units_v3`'s
  unit definitions were built using the TESTIS reads. §4l's 87% edge-recovery figure is about O1's
  family *topology* transporting across substrates, not a guarantee that every unit's exact boundary
  realigns equally well to fibroblast reads for O3's specific consistent-mismatch-site statistic. Before
  trusting any differential candidate, check that the overwhelming majority of tested pairs land in
  "no signal" (both `None`) — if the split looks anywhere near even, the catalog itself is the confound,
  not individual biology.
- Depth-matching (`fibro_ds.bam` already downsampled to `GGO_ds.bam`'s depth) controls for library-size
  differences but not for tissue-specific library complexity/duplication-rate differences, which are not
  independently checked here.

## Open questions for the plan (not blocking this design, but worth flagging)

- Exact wall-clock cost of the added realignment step on a genome-wide sweep is unmeasured; the Python
  script already has a `--budget` (seconds) escape valve the Rust port may or may not need — decide once
  a real timing number exists.
- Whether `lgamma` in `allele_specific_junctions.rs` should be made `pub(crate)` and shared, or
  duplicated — a two-line decision, deferred to the implementation task rather than blocking this design.
