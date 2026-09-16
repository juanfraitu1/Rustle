# Read-Dedup Fix + `--gtf-refine` Bundle, Validated on Held-Out chr17 — Design

**Status**: DESIGNED, not implemented. Written 2026-09-16 from the chr20 error-pattern dissection
(`/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/diagnostics/`, 4-agent workflow, headline counts
independently recounted). chr20 is the DEVELOPMENT substrate for everything below; chr17 is held out.

## Why

On chr20 our `copy_assign --gtf` assembly is in StringTie/FLAIR's range but loses on precision (transcript
Pr 35.6 vs StringTie 47.1). The dissection found the losses are NOT junction errors (multi-exon 'j' 287 vs
StringTie 285; NNC 166 vs 177) — every read-consensus / splice-site-snapping variant was net negative
because nearby junction variants are 25/26 canonical real alternative sites. Instead:

- **Precision** is lost to three missing, annotation-free collapse/filter steps that FLAIR, StringTie and
  isoseq collapse all have: single-exon strand is a hard-coded `'+'` placeholder (all 203/203 mono models);
  own-subset (truncated-fragment) models are never removed (94 multi-exon 'c' vs StringTie 32); single-exon
  models inside our own spliced exons / at spliced-read-dominated loci are never removed.
- **Sensitivity** is lost to support COUNTING, not junction accuracy: 48/56 of the reference transcripts that
  StringTie or FLAIR recover and we miss have exactly one full-length read; 28 have compatible truncated reads
  that the other tools count and we do not (27 of those 28 have exactly one full-length read).
- **A real bug** (`src/bin/copy_assign.rs` ~2403-2410): the primary-read dedup meant to drop a read returned
  by two adjacent windows keys on placement `(chrom, ref_start, ref_end, introns)` and runs WITHIN each window,
  so distinct molecules with identical placement collapse to one. chr20 (one window): 25,341 primaries →
  11,755 (53.6% of real molecules dropped; only 474 are same-UMI duplicates). Its own comment states the
  intended key is per-molecule.

## Scope

1. Fix the dedup bug, with an escape hatch, and MEASURE its effect on O2 (`--families`) and de novo
   `copy_assign` before merge.
2. Add `--gtf-refine`, an opt-in bundle of four annotation-free components for the `--gtf` path only.
3. Pre-register chr17 decision rules, reproduce the chr20 simulation numbers exactly with the Rust
   implementation (fidelity), then run the held-out chr17 evaluation (gffcompare + SQANTI3, with StringTie and
   FLAIR on chr17 for context).

**Non-goals**: splice-site correction of any kind (measured net negative on chr20 — rejected); raising or
lowering support floors (3-read floor lost 41 '='; singleton floor dropped transcript Pr to 18.6 — rejected);
`RUSTLE_JUNCTION_MAJORITY` (documented unmeasured chr16 multi-copy harm); any change to
`pass1_skeletons_robust`, `detect_and_assign`, or `gw_family_catalog`.

## Part 1 — dedup fix

**Change** (`src/bin/copy_assign.rs`, the `seen.insert(...)` loop over `wins`): dedup a primary read only
against keys kept from EARLIER windows, never against reads of the current window. After each window, add its
kept keys to the seen-set.

Why this is exact without read names (PrimaryRead has no name; adding one would cost memory genome-wide): a
cross-window key match implies identical `(chrom, start, end, chain)`, hence the read overlaps the earlier
window too and was already returned by it — so cross-window matches are exactly the boundary duplicates, and
two distinct molecules with identical placement are both kept from the first window that returns them.

**Escape hatch**: `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` restores the old within-window placement dedup
byte-for-byte (precedent: `RUSTLE_EDGE_CORE=poa`), so every past ledger number stays reproducible.

**Default: fixed.** ⚠ DECISION POINT FOR THE USER at merge time (not before): the measurement below is shown
first; if it reveals changes the user does not want in the default, flipping the default is one line.

**Pre-registered impact predictions** (checked before merge):
- **O2 (`--families`)**: `assignments.tsv`, `quant.tsv`, `families.tsv`, `family_join.tsv` BYTE-IDENTICAL
  fixed vs legacy on the real `mec/run_psv.sh` substrate (`mec/batch.copies.tsv`, `mec/regions.txt`,
  `npip_cat/npip3.bam`, `npip_cat/arm_f2/cat.copies.fa`). Reason: under `--families` the skeleton front end is
  empty and rescue is skipped; the reads used for assignment (`bam_reads`) are deduplicated separately and
  correctly, by `(name, chrom, ref_start)`.
  **If not byte-identical: STOP and report before any further task.**
- **De novo `copy_assign` (no `--families`)** on the same regions: expected to CHANGE (pass-1 floors, mischain
  junction support, thin-loci rescue all see more reads). Report families/copies/assignment-status deltas;
  descriptive, no pass/fail.
- **chr20 `--gtf`** (no refine): exact reproduction of simulation arm A1 — see Fidelity.

## Part 2 — `--gtf-refine`

CLI: `--gtf-refine <list>` (comma-separated: `strand`, `subset`, `mono`, `fragsupport`, or `all`; default
empty). Affects ONLY the `if args.gtf { ... }` block in `src/bin/copy_assign.rs` (the site that calls
`pass1_skeletons(&primary, ...)`). Empty list ⇒ byte-identical output. No component reads the annotation, so
the bundle is usable in both de novo and guided modes. Pure logic lives in a new lib module
`src/rustle/vg_family/gtf_refine.rs` (registered OPT-IN in `docs/MODULE_STATUS.md`), unit-tested with
synthetic data.

Pipeline order when enabled: (1) pass-1 with `fragsupport` → (2) `assemble_gate_with(..., use_read_strand =
strand enabled, strand_margin = 0.90)` → (3) post-filters `subset`, then `mono` → (4) existing
`collapse_loci_groups` and GTF writing. Thresholds below are FIXED (chosen on chr20; not to be tuned on chr17).

### `fragsupport` — 3'-anchored fragment support (chr20 arm B3)

Spliced primary reads grouped by exact `(chrom, intron chain)` as today. Candidate chains = chains with ≥1
exact read. A spliced read R with chain r (≥1 intron; unspliced reads never count) is a FRAGMENT of candidate
F iff all hold:
- `len(F) > len(r)` and `F[i:i+m] == r` for some i (contiguous sub-chain);
- intron compatibility: if `i > 0`, `R.ref_start >= F[i-1].acceptor`; if `i+m < len(F)`,
  `R.ref_end <= F[i+m].donor`;
- 3'-anchored on F's motif strand (strict all-canonical chain strand; F with no strict strand gets no
  fragments): `'+'` requires `i+m == len(F)`; `'-'` requires `i == 0`;
- candidates are looked up via chains containing r's first intron.
Assign-or-abstain: R counts only if it is a fragment of EXACTLY ONE candidate (no 1/k splitting). Emit F if
`exact(F) + fragments(F) >= 2`; `start`/`end` from exact reads only (min start / max end); `n_reads` =
`exact + fragments`. Unspliced clustering unchanged.

### `strand` — read-orientation strand for single-exon models

Use the existing machinery via `assemble_gate_with(..., true, 0.90)` for the `--gtf` call only (vote over the
model's own pass-1 cluster reads; `'-'` iff a strict majority reverse and majority fraction ≥ 0.90, else
`'+'` placeholder). No env var involved; no other `assemble_gate` call changes.

### `subset` — own-subset removal (chr20 `sub_own_ovh5`)

Over the multi-exon models after the gate, computed against the whole pre-filter set at once: remove T if
some other model S on the same chrom and strand has `S.introns[p:p+k] == T.introns` (k = len(T.introns) ≥ 1)
and both ends fit: left overhang `lo = S.exon[p].start - T.start` satisfies `lo <= 0 || (p != 0 && lo <= 5)`,
right overhang `ro = T.end - S.exon[p+k].end` satisfies `ro <= 0 || (p+k != len(S.introns) && ro <= 5)`.
No support condition.

### `mono` — single-exon removal (chr20 P1 mono-in-exon + P5 spliced dominance)

Applied after `subset`, to single-exon models only:
1. **Inside own spliced exon (strand-aware)**: remove M if some remaining model with ≥2 exons, same strand as
   M (M's strand after `strand`, if enabled), has an exon `[s, e)` with `s <= M.start && e >= M.end`.
2. **Spliced-read dominated**: over ALL primary reads gathered for the region (post-dedup-fix), with
   `U` = unspliced primaries overlapping M, `SI` = spliced primaries with an intron covering ≥ 50% of M's
   length, `SE` = spliced primaries on M's strand (read strand = FLAG 0x10) with an aligned exon block
   covering ≥ 50% of M's length: remove M if `SI >= U || SE >= U`.
   (Known deviation from the chr20 simulation, which set spliced-read strand as FLAG 0x10 flipped by
   minimap2's `ts` tag; `PrimaryRead` carries no `ts`. FLAG agrees with motif strand in 99.95% of chr20
   spliced primaries, so the effect is expected to be ~0 and is checked in Fidelity.)

## Part 3 — fidelity on chr20, then held-out chr17

### Fidelity anchors (chr20, Rust implementation; gffcompare vs `chr20_ref.gtf`)

EXACT (every deviation must be explained down to individual transcripts before any chr17 file is generated;
if the explanation is an error in the Python simulation rather than the Rust code, record that explicitly —
do not bend the implementation toward a wrong simulation):
| run | queries | '=' | transcript Sn/Pr | intron-chain Sn/Pr |
|---|---|---|---|---|
| legacy dedup, no refine | 976 | 347 | 7.6 / 35.6 | 8.0 / 44.6 |
| fixed dedup, no refine (A1) | 1022 | 350 | 7.7 / 34.2 | 8.1 / 42.9 |
| fixed dedup + `fragsupport` (B3) | 1054 | 362 | 7.9 / 34.3 | 8.4 / 42.7 |
| legacy dedup + `subset` | 913 | 347 | 7.6 / 38.0 | 8.0 / 48.6 |

(A1, the legacy baseline and legacy + `subset` were independently recomputed during the dissection; B3 rests
on the validated Python replica only.)

APPROXIMATE (deviations allowed only if explained by the named differences: vote over deduplicated cluster
reads under legacy dedup vs all overlapping unspliced primaries in the simulation; FLAG-only spliced-read
strand):
| run | queries | '=' | notes |
|---|---|---|---|
| legacy dedup + `strand` | 976 | 349 | locus Pr 50.4 |
| legacy dedup + `strand,subset,mono` (P5) | 820 | 349 | transcript Pr 42.6, locus Pr 52.2 |

Also reported (development numbers, not validation): fixed dedup + `all` on chr20 — never simulated as a
combination.

### Held-out chr17

Substrate: `human_testis.t2t.bam` restricted to chr17 (84,276,897 bp; 462,463 records), CHM13 genome, chr17
subset of `chm13v2.0_RefSeq_full.gff.gz` (same prep as `bench/prep_chr20_ref.sh`, incl. its GFF3 resort fix).
Chosen over chrX (twice the size; testis-specific multi-copy cancer-testis families would confound an
ordinary-assembly test). Note: chr17 carries the TBC1D3 / 17q21.31 duplications this project has studied as
FAMILIES — not for assembly tuning, so it remains held out for this purpose. Outputs under
`/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17/`.

Runs: ours BASELINE (fixed dedup, no refine) and ours BUNDLE (fixed dedup, `--gtf-refine all`) — both
gffcompare + SQANTI3; ours legacy (current shipped behavior, context); ablations (each component alone,
gffcompare only); StringTie `-L` and FLAIR (unguided, same workarounds as chr20) — both gffcompare + SQANTI3,
context only.

**Pre-registered decision rules** (bundle vs baseline; must be committed in
`docs/PREREG_gtf_refine_chr17_2026-09-16.md` BEFORE any chr17 file is generated):
- **E1 matches**: distinct reference transcripts with '=' — bundle ≥ baseline.
- **E2 precision**: transcript-level Pr AND intron-chain-level Pr — bundle > baseline (strict).
- **E3 collateral**: references '=' in baseline but not in bundle ≤ 1% of baseline's '=' count.
- **E4 novel junctions**: distinct canonical novel junctions (SQANTI3 `junction_category == novel`, canonical,
  any transcript) — bundle ≥ baseline. By construction no component removes a junction not retained in a kept
  model; a violation indicates an implementation bug.
- **Verdict**: SUPPORTED if E1–E4 all hold. PARTIAL if E3 and E4 hold and exactly one of E1/E2 holds.
  REFUTED otherwise. Ablations, SQANTI3 category tables, the multicopy-tagged subset, and StringTie/FLAIR are
  descriptive only.
No threshold, window, or component may be changed after any chr17 number is seen; a changed rule is a new,
separately pre-registered experiment on a different substrate.

## Testing

- Unit tests (`gtf_refine.rs`): fragment compatibility (suffix vs prefix by strand, intron-compatibility
  rejection at both ends, ambiguous fragment abstains, unspliced reads ignored, emission threshold), subset
  overhang rules (terminal vs non-terminal container exon, 5 bp boundary, strand mismatch never removes),
  mono containment (strand-aware), spliced-dominance counts.
- Unit test for the dedup rule: two distinct reads with identical placement in one window are both kept; the
  same read returned by two windows is kept once.
- Integration (`tests/copy_assign_families.rs` conventions): `--gtf-refine` empty ⇒ byte-identical outputs;
  `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` ⇒ byte-identical to the pre-fix binary (compare against output generated
  from the base commit). The fixture's `--gtf` emits 0 rows (known), so ON-state behavior is proven by the
  unit tests and the chr20 exact fidelity anchors, and this limitation is stated in the test's doc comment.

## Global constraints

- Default outputs change ONLY through the dedup fix; `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` reproduces the pre-fix
  binary byte-for-byte on every output.
- `--gtf-refine` empty ⇒ byte-identical to the fixed-dedup binary.
- Never modify `pass1_skeletons_robust`, `detect_and_assign`, or `gw_family_catalog`; `--gtf-refine` touches
  only the `if args.gtf` block.
- No component uses the annotation. Thresholds exactly as written above.
- O2 byte-identity is a hard stop condition.
- WSL2 rules: builds with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target`, one heavy run at a
  time in the foreground, output redirected to files, never `pkill -f`, big outputs under `/mnt/linuxdisk`.
