# PREREG — held-out chr17 evaluation of `--gtf-refine` (2026-09-16, written before any chr17 file exists)

Step 1 of this task confirmed `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17` does not exist. No chr17
BAM, FASTA, GTF, gffcompare output, or SQANTI3 output has been generated, subsetted, or looked at. This
commit is the timestamp that makes the chr17 test blind: **no threshold, component, or rule below may change
after any chr17 number is seen.**

## Question

On held-out chr17 (never used to develop the dedup fix or any `--gtf-refine` component), does the
`--gtf-refine all` bundle (fixed dedup + `strand` + `subset` + `mono` + `fragsupport` + `tss`) improve our
de novo `copy_assign --gtf` assembly over the fixed-dedup baseline, without the precision/junction losses
that `--gtf` currently pays relative to StringTie and FLAIR — and does the `tss` component's 5'-end rule move
transcript starts closer to the true TSS?

## Substrate

- Reads: `human_testis.t2t.bam` restricted to chr17 (84,276,897 bp; 462,463 records).
- Genome: CHM13 (`chm13v2.0.fa`).
- Annotation: chr17 subset of `chm13v2.0_RefSeq_full.gff.gz`, prepared exactly as
  `bench/prep_chr20_ref.sh` prepared the chr20 reference (chr17-filtered GFF3, then the same GFF3→GTF
  conversion including its GFF3 resort fix — the upstream RefSeq GFF3 lists a feature's leftmost exon before
  its own gene/transcript record whenever they share a start coordinate, which silently collapses `gene_id`
  to `transcript_id` in gffread's single-pass parser unless the file is stably re-sorted by
  `(chrom, start, feature-rank)` first; the same bug and the same fix apply verbatim to chr17).
- Outputs live under `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17/`.
- chr17 was chosen over chrX (twice the size; testis-specific multi-copy cancer-testis gene families would
  confound an ordinary-assembly test). chr17 carries the TBC1D3 / 17q21.31 duplications this project has
  separately studied as gene FAMILIES; that prior study is not used to tune assembly here, and chr17 is held
  out for this assembly-evaluation purpose only.

## Arms

- **BASELINE**: fixed dedup (Part 1 of the spec), no `--gtf-refine` (empty list).
- **BUNDLE**: fixed dedup, `--gtf-refine all` (all five components: `strand`, `subset`, `mono`,
  `fragsupport`, `tss`).
- **ABL_TSS**: fixed dedup, `--gtf-refine tss` only — scored for E5 only, against BASELINE, so that only the
  5' end differs between the two classification files compared.
- Descriptive runs (not scored against E1–E5, context only): **legacy** (current shipped behavior: legacy
  dedup, no refine); the four other single-component ablations (`strand` alone, `subset` alone, `mono`
  alone, `fragsupport` alone — gffcompare only); **StringTie `-L`** and **FLAIR** (unguided, same workarounds
  used on chr20) — both gffcompare + SQANTI3, context only.
- BASELINE and BUNDLE are scored with gffcompare + SQANTI3. Ablations other than `abl_tss` get gffcompare
  only. `abl_tss` additionally gets SQANTI3 because E5 is scored on it.

## Frozen component rules (verbatim from `docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md`, Part 2)

CLI: `--gtf-refine <list>` (comma-separated: `strand`, `subset`, `mono`, `fragsupport`, `tss`, or `all`;
default empty; `all` = all five — `tss` added by the 2026-09-16 addendum below). Affects ONLY the
`if args.gtf { ... }` block in `src/bin/copy_assign.rs` (the site that calls
`pass1_skeletons(&primary, ...)`). Empty list ⇒ byte-identical output. No component reads the annotation, so
the bundle is usable in both de novo and guided modes. Pure logic lives in a new lib module
`src/rustle/vg_family/gtf_refine.rs` (registered OPT-IN in `docs/MODULE_STATUS.md`), unit-tested with
synthetic data.

Pipeline order when enabled: (1) pass-1 with `fragsupport` → (2) `assemble_gate_with(..., use_read_strand =
strand enabled, strand_margin = 0.90)` → (3) post-filters `subset`, then `mono` → (3b) `tss` 5'-end
refinement → (4) existing `collapse_loci_groups` and GTF writing. Thresholds below are FIXED (chosen on
chr20; not to be tuned on chr17).

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

### `tss` — densest 5'-start window (ADDENDUM 2026-09-16, user-approved before any chr17 data); W = 50 (frozen)

Why: the transcript 5' end is currently the most extreme exact read, so one outlier read sets it. On chr20
(legacy run, SQANTI3 full-splice matches, multi-exon) 80/345 models start >50 bp upstream of the matched
reference TSS vs 7/269 for FLAIR, and the overshoot grows with support (≥10 reads: 53 overshoot vs 70 within
50 bp; 2 reads: 5 vs 58) — the signature of an extreme-read boundary. gffcompare transcript/intron-chain
metrics ignore multi-exon ends, so E1–E4 cannot reward this; it gets its own endpoint (E5).

Rule, applied after `mono` (step 3b), to each model with ≥1 intron and strand `'+'` or `'-'` (single-exon
models and all 3' ends unchanged):
- Exact reads = primary reads (the same post-dedup-fix `primary` set) on the model's chrom whose intron chain
  equals the model's exactly. Fragments never count (their 5' ends are truncated by definition). A model with
  0 exact reads is left unchanged.
- 5' ends: `'+'` → `ref_start`; `'-'` → `ref_end`. Sort ascending as `v`.
- `'+'`: for each i, `c_i = #{x in v : v[i] <= x <= v[i] + W}`; choose the maximal `c_i`, ties → smallest i
  (most upstream); new `start = v[i]`.
- `'-'`: for each j, `c_j = #{x in v : v[j] - W <= x <= v[j]}`; choose the maximal `c_j`, ties → largest j
  (most upstream on `'-'`); new `end = v[j]`.
- Consequences (not extra rules): the 5' end can only move downstream of the most extreme read or stay; it
  stays strictly upstream of the first donor (`'+'`) / downstream of the last acceptor (`'-'`) because every
  exact read does; a model with 1–2 exact reads never changes (two reads within W: the upstream window holds
  both; more than W apart: a 1–1 tie keeps the most upstream read), so small models are untouched without any
  read-count cutoff. A change needs ≥2 reads clustered downstream of ≥1 outlier.

`W` was chosen on chr20 (the development substrate), never on chr17, and is now frozen as a Rust constant:
`pub const TSS_WINDOW_BP: u64 = 50;` (`src/rustle/vg_family/gtf_refine.rs:240`). Selection record
(`bench/CHR20_ASSEMBLER_COMPARISON.md`, "Follow-up: `tss` window selection" section, Task 7): a Python
simulation over W ∈ {10, 25, 50, 100} on the fixed-dedup, no-refine chr20 model set (A1), selecting the W
maximizing the count of multi-exon models (with a chain matching ≥1 chr20 reference transcript) whose new 5'
end is within 50 bp of that reference's TSS, ties → smaller W. Result: `n_within50` = 246 (none) / 297 (W10) /
296 (W25) / **303 (W50)** / 278 (W100) of 348 chain-matched models — W = 50 chosen. The Rust `tss` component
was then checked (`bench/CHR20_ASSEMBLER_COMPARISON.md`, "Follow-up: chr20 fidelity + SQANTI3 check for
`tss`" section, Task 9) to reproduce the chosen-W simulation's GTF exactly, transcript-by-transcript
(1022/1022 match, 0 only-in-rust, 0 only-in-sim).

## Development record (chr20; cited, not re-derived)

- All Rust fidelity anchors (EXACT and APPROXIMATE, dedup fix + `strand`/`subset`/`mono`/`fragsupport`)
  reproduced the Python simulation exactly on real chr20 data, including per-transcript checks (0 deviating
  transcripts in every arm) — `bench/CHR20_ASSEMBLER_COMPARISON.md`, "Follow-up: Rust fidelity of the dedup
  fix and `--gtf-refine`" section (Task 6/8).
- `fixed_tss` (fixed dedup + `--gtf-refine tss`) reproduced its Python simulation transcript-by-transcript:
  1022/1022 transcripts identical (chrom, strand, intron chain, first-exon start, last-exon end); 0 only in
  Rust, 0 only in simulation — `bench/CHR20_ASSEMBLER_COMPARISON.md`, Task 9 section.
- SQANTI3 on chr20, `fixed_none` vs `fixed_tss` (full-splice-match multi-exon models, n = 348): `p` (fraction
  `|diff_to_TSS| <= 50`) 0.7069 → 0.8707; `g` (count `|diff_to_gene_TSS| <= 50`) 268 → 319 — same section,
  descriptive (chr20 is development, not the pre-registered test).

## Descriptive vs decision-bearing

Only BASELINE vs BUNDLE (E1–E4) and BASELINE vs ABL_TSS (E5) carry pre-registered pass/fail verdicts.
Everything else on chr17 — legacy, the four single-component ablations, SQANTI3 category tables, the
multicopy-tagged subset, and StringTie/FLAIR — is descriptive context only, reported alongside the verdicts,
never substituted for them.

## Scoring conventions

- **`'='` match**: a distinct `ref_id` with `class_code == '='` in the gffcompare `.tmap` output.
- **Precision**: the one-decimal-precision value read from the gffcompare `.stats` line (`Transcript level:`
  or `Intron chain level:`, second `|`-delimited number). A tie at one decimal place is a FAILURE of a
  strict (`>`) comparison — e.g. E2 requires `bundle_precision > baseline_precision` at that one-decimal
  value; equal one-decimal values fail E2, even if the underlying unrounded precision differs.
- **Novel canonical junction**: a distinct `(chrom, strand, genomic_start_coord, genomic_end_coord)` tuple
  with `junction_category == novel` and `canonical == canonical` in SQANTI3's `_junctions.txt`, pooled over
  all transcripts in the run.
- **E5 row filter**: SQANTI3 `_classification.txt` rows with `structural_category == full-splice_match` AND
  `subcategory != mono-exon`. `p` is the UNROUNDED fraction `within / n` (not the one-decimal precision
  convention used for E1–E4). A row whose `diff_to_TSS` or `diff_to_gene_TSS` field is non-numeric counts
  toward `n` but not toward `within`/`g` for that field (i.e. it is excluded only from the within/gene
  numerators, not from the denominator).
- The verdict is computed by `bench/gtf_refine_verdict.py` (this commit) — no other script or manual
  arithmetic decides E1–E5.

## Decision rules (frozen)

Verbatim from the spec's Part 3 "Held-out chr17" section (bundle vs baseline; the numbers below are computed
by `bench/gtf_refine_verdict.py`'s `verdict()`, this commit):
- **E1 matches**: distinct reference transcripts with '=' — bundle ≥ baseline.
- **E2 precision**: transcript-level Pr AND intron-chain-level Pr — bundle > baseline (strict).
- **E3 collateral**: references '=' in baseline but not in bundle ≤ 1% of baseline's '=' count.
- **E4 novel junctions**: distinct canonical novel junctions (SQANTI3 `junction_category == novel`, canonical,
  any transcript) — bundle ≥ baseline. By construction no component removes a junction not retained in a kept
  model; a violation indicates an implementation bug.
- **Verdict**: SUPPORTED if E1–E4 all hold. PARTIAL if E3 and E4 hold and exactly one of E1/E2 holds.
  REFUTED otherwise. Ablations, SQANTI3 category tables, the multicopy-tagged subset, and StringTie/FLAIR are
  descriptive only.
- **E5 5' ends (addendum; a separate claim with its own verdict, scored `abl_tss` vs BASELINE so only 5'
  ends differ)**: over SQANTI3 rows with `structural_category == full-splice_match` and `subcategory !=
  mono-exon`, let `p` = fraction with `|diff_to_TSS| <= 50` (unrounded) and `g` = count with
  `|diff_to_gene_TSS| <= 50`. E5 holds iff `p(abl_tss) > p(baseline)` AND `g(abl_tss) >= g(baseline)` (the
  guard keeps real alternative TSSs from being penalised). **TSS verdict**: SUPPORTED if E5 holds, else
  REFUTED. Reported next to, not merged into, the E1–E4 verdict. 50 bp is SQANTI3's own reference_match
  tolerance, not a new threshold. (The E5/TSS-verdict numbers above are computed by
  `bench/gtf_refine_verdict.py`'s `tss_verdict()`, this commit.)

No threshold, component, or rule may change after any chr17 number is seen; a changed rule is a new,
separately pre-registered experiment on a different substrate.

Note on E2 for the self-review record: "strict" applies independently to BOTH the transcript-level and
intron-chain-level precision comparisons (`bun_pr[0] > base_pr[0] and bun_pr[1] > base_pr[1]` in
`gtf_refine_verdict.py`) — a tie or a loss on either one fails E2, even if the other precision improves.

## Known limits (declared now)

- `tss` can only move a 5' end downstream (toward the read cluster) or leave it unchanged, never upstream of
  the most extreme exact read. On chr20, 22 of 348 chain-matched models had the reference TSS upstream of
  every exact read, so `tss` cannot reach the true TSS for them regardless of W (ceiling 326/348, not
  348/348) — `progress.md` line 72 (Task 7 observation).
- `tss` can move a well-supported model AWAY from its true TSS when most reads are 5'-truncated: chr20
  example `DN_chr20_37420188_3` (132 reads) moved from 8 bp to 168 bp from the reference TSS — the majority
  of 5'-truncated reads outvoted the full-length minority. This is 1 of 6 `moved_out` cases at W = 50 on
  chr20 (`bench/CHR20_ASSEMBLER_COMPARISON.md`, Task 7 table; `progress.md` line 72).
- E1–E4 are scored against a single annotation (RefSeq/CHM13). A component that assembles a real but
  unannotated isoform correctly is scored as a miss, not a match; the ground-truth ceiling this implies is
  not measured here.
- chr17 contains the TBC1D3 / 17q21.31 duplications this project has separately studied as gene FAMILIES
  (see the project's family-definition work). That study is not used to tune this assembly evaluation, and
  no family-level result from it is imported into the E1–E5 verdicts.
- `mono`'s spliced-read strand uses FLAG only (0x10), not FLAG⊕`ts` as the chr20 simulation did; this is
  checked to have ~0 effect on chr20 (99.95% FLAG/motif-strand agreement) but is not re-derived on chr17
  before scoring — a real, larger effect on chr17 would show up as an unexplained E1–E4 deviation, not a
  changed rule.
- `fixed_all` (all five components together) was never itself simulated in Python as a combination on chr20;
  only the individual components and, for `tss`, one additional pairwise arm (`fixed_tss`) were validated
  against a simulation. BUNDLE on chr17 therefore carries less direct fidelity evidence than any single
  component does.
