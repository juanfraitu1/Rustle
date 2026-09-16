# PREREG — held-out chr17 evaluation of `--gtf-refine` (2026-09-16, written before any chr17 file exists)

Confirmed with `ls -d /mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17` that the directory does not exist.
No chr17 BAM, FASTA, GTF, gffcompare output, or SQANTI3 output has been generated, subsetted, or looked at.
This commit is the timestamp that makes the chr17 test blind: **no threshold, window, or component below may
change after any chr17 number is seen.**

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
  dedup, no refine) — gffcompare only; the four other single-component ablations (`strand` alone, `subset`
  alone, `mono` alone, `fragsupport` alone) — gffcompare only; **StringTie `-L`** and **FLAIR** (unguided,
  same workarounds used on chr20) — both gffcompare + SQANTI3, context only.
- BASELINE and BUNDLE are scored with gffcompare + SQANTI3. `legacy` and the four other ablations
  (`abl_strand`, `abl_subset`, `abl_mono`, `abl_fragsupport`) get gffcompare only. `abl_tss` additionally
  gets SQANTI3 because E5 is scored on it.
- StringTie and FLAIR are context only and never enter any E1–E5 or TSS verdict computation, at any point.

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

Addition not in the spec's Part 2 prose above (found in `copy_assign.rs:2716-2732` and named here for
completeness, since it sits inside the same `if args.gtf` block between step (1) and step (2)): an opt-in
`RUSTLE_JUNCTION_FUZZ_BP` fuzzy intron-chain merge (default off/unset = no-op, `merge_fuzzy_skeletons` with
tolerance 0). This step is **OFF (unset) for every arm of this test** — it was pre-registered and measured
net-negative on chr20 at the 672 bp tolerance (`bench/CHR20_ASSEMBLER_COMPARISON.md`, "Follow-up: fuzzy
junction merging at the pre-registered tolerance (672bp) — NEGATIVE" section) and is not one of the five
`--gtf-refine` components. See "Run protocol (frozen)" below for how it (and every other un-named `RUSTLE_*`
variable) is forced off.

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

Implementation detail confirmed against `gtf_refine.rs`'s `subset_removals` (not in the spec prose, added here
for precision): S must have strictly more introns than T (`s.introns.len() > k`, i.e. `len(S.introns) >
len(T.introns)`) — two models with identical intron chains never remove each other under this rule.

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
simulation over W ∈ {10, 25, 50, 100} on the fixed-dedup, no-refine chr20 model set (A1). Selection metric:
among A1 multi-exon models whose (strand, intron chain) equals ≥1 chr20 reference transcript's, count models
with `min |model 5' end − ref TSS| <= 50` over those references; choose the W maximizing the count, ties →
smaller W. Result: `n_within50` = 246 (none) / 297 (W10) / 296 (W25) / **303 (W50)** / 278 (W100) of 348
chain-matched models — W = 50 chosen. The Rust `tss` component was then checked (`bench/CHR20_ASSEMBLER_COMPARISON.md`,
"Follow-up: chr20 fidelity + SQANTI3 check for `tss`" section, Task 9) to reproduce the chosen-W simulation's
GTF exactly, transcript-by-transcript (1022/1022 match, 0 only-in-rust, 0 only-in-sim).

## Development record (chr20; cited, not re-derived)

- All Rust fidelity anchors (EXACT and APPROXIMATE, dedup fix + `strand`/`subset`/`mono`/`fragsupport`)
  reproduced the Python simulation exactly on real chr20 data, including per-transcript checks (0 deviating
  transcripts in every arm, including the 3-component `legacy_strand_subset_mono` (P5) arm at 820/820) —
  `bench/CHR20_ASSEMBLER_COMPARISON.md`, "Follow-up: Rust fidelity of the dedup fix and `--gtf-refine`"
  section (Task 6). Those per-transcript checks were keyed on `(chrom, strand, intron chain)` for multi-exon
  models (start/end were not yet part of any rule under test); `fixed_tss`'s check (next bullet) added the
  5' coordinate to the key because `tss` is the first component that can change it.
- `fixed_tss` (fixed dedup + `--gtf-refine tss`) reproduced its Python simulation transcript-by-transcript,
  keyed on `(chrom, strand, intron chain, first-exon start, last-exon end)`: 1022/1022 transcripts identical;
  0 only in Rust, 0 only in simulation — `bench/CHR20_ASSEMBLER_COMPARISON.md`, Task 9 section.
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
- **E5 zero-row edge case**: if `n = 0` for either the BASELINE or the ABL_TSS arm (no full-splice-match
  multi-exon SQANTI3 rows at all), `p = 0.0` for that arm (`tss_metrics()`'s `within / n if n else 0.0`) and
  the TSS verdict is REFUTED. This is declared now, before any chr17 number exists, so that a degenerate
  denominator is never argued into a SUPPORTED reading after the fact.

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

No threshold, window, or component may be changed after any chr17 number is seen; a changed rule is a new,
separately pre-registered experiment on a different substrate. (Spec's exact wording, Part 3.)

Note on E2 for the self-review record: "strict" applies independently to BOTH the transcript-level and
intron-chain-level precision comparisons (`bun_pr[0] > base_pr[0] and bun_pr[1] > base_pr[1]` in
`gtf_refine_verdict.py`) — a tie or a loss on either one fails E2, even if the other precision improves.

## Run protocol (frozen)

This commit freezes the decision rules; it does not by itself freeze how the runs are executed and scored.
This section closes that gap. The five scripts below are committed together with this PREREG (same commit),
so their content is frozen at the same instant as the rules above.

**Scripts** (all under `bench/`, chromosome-parameterized siblings of the chr20 bakeoff scripts — the chr20
scripts themselves are never modified):
- `bench/prep_chrom_ref.sh` — substrate prep (BAM/FASTA/reference-GTF extraction, same GFF3 resort fix as
  `bench/prep_chr20_ref.sh`).
- `bench/bakeoff_chrom_ours.sh` — runs `copy_assign --gtf` for one arm/label.
- `bench/bakeoff_chrom_stringtie.sh` — StringTie `-L`, unguided, context only.
- `bench/bakeoff_chrom_flair.sh` — FLAIR, unguided (same chr20 PATH/`flair correct`-skip workarounds),
  context only.
- `bench/chrom_score.sh` — gffcompare (+ optional SQANTI3) for one label.

**Exact commands, chr17** (`CHROM=chr17`; `W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17`):
```
bash bench/prep_chrom_ref.sh chr17

bash bench/bakeoff_chrom_ours.sh chr17 baseline
bash bench/bakeoff_chrom_ours.sh chr17 bundle --gtf-refine all
bash bench/bakeoff_chrom_ours.sh chr17 abl_tss --gtf-refine tss
RUSTLE_LEGACY_PLACEMENT_DEDUP=1 bash bench/bakeoff_chrom_ours.sh chr17 legacy
bash bench/bakeoff_chrom_ours.sh chr17 abl_strand --gtf-refine strand
bash bench/bakeoff_chrom_ours.sh chr17 abl_subset --gtf-refine subset
bash bench/bakeoff_chrom_ours.sh chr17 abl_mono --gtf-refine mono
bash bench/bakeoff_chrom_ours.sh chr17 abl_fragsupport --gtf-refine fragsupport

bash bench/bakeoff_chrom_stringtie.sh chr17
bash bench/bakeoff_chrom_flair.sh chr17

bash bench/chrom_score.sh chr17 baseline          $W/baseline/ours.gtf         --sqanti
bash bench/chrom_score.sh chr17 bundle            $W/bundle/ours.gtf           --sqanti
bash bench/chrom_score.sh chr17 abl_tss           $W/abl_tss/ours.gtf          --sqanti
bash bench/chrom_score.sh chr17 legacy            $W/legacy/ours.gtf
bash bench/chrom_score.sh chr17 abl_strand        $W/abl_strand/ours.gtf
bash bench/chrom_score.sh chr17 abl_subset        $W/abl_subset/ours.gtf
bash bench/chrom_score.sh chr17 abl_mono          $W/abl_mono/ours.gtf
bash bench/chrom_score.sh chr17 abl_fragsupport   $W/abl_fragsupport/ours.gtf
bash bench/chrom_score.sh chr17 stringtie         $W/stringtie/st.gtf          --sqanti
bash bench/chrom_score.sh chr17 flair             $W/flair/flair.isoforms.gtf  --sqanti
```
BASELINE: no `--gtf-refine` flag, `RUSTLE_LEGACY_PLACEMENT_DEDUP` unset (fixed dedup, the Part 1 default).
BUNDLE: `--gtf-refine all`. ABL_TSS: `--gtf-refine tss` only, scored for E5 only. legacy: fixed dedup binary
run with `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` (the ONLY named exception to the blanket `RUSTLE_*` unset below).
The four other ablations run one named component each. SQANTI3 (`--sqanti`) is run for exactly five labels:
`baseline`, `bundle`, `abl_tss`, `stringtie`, `flair` — matching the spec's Part 3 Runs paragraph and the
plan's Task 11 file list. `legacy`, `abl_strand`, `abl_subset`, `abl_mono`, `abl_fragsupport` get gffcompare
only.

**Environment policy.** `bench/bakeoff_chrom_ours.sh` force-unsets every environment variable whose name
starts with `RUSTLE_` except `RUSTLE_LEGACY_PLACEMENT_DEDUP`, unconditionally, before running the binary —
regardless of what the calling shell happens to have exported. This is what keeps `RUSTLE_TSS_SNAP` (and its
`_WINDOW`/`_FRAC` companions), `RUSTLE_FOOTPRINT_NODES`, `RUSTLE_JUNCTION_FUZZ_BP`,
`RUSTLE_JUNCTION_MAJORITY`, `RUSTLE_JUNCTION_NC_MAX_BP`, and any future `RUSTLE_*` variable off for every
scored arm without having to enumerate them individually in this document.

**Binary provenance.** The binary must be built from a tree whose `src/`, `Cargo.toml`, and `Cargo.lock` are
byte-identical to commit `3b096012` (the commit that carries this PREREG and `bench/gtf_refine_verdict.py`)
— chosen as the base because it is the commit that froze the rules being tested; a later commit is
acceptable ONLY if it touches none of `src/`, `Cargo.toml`, `Cargo.lock`. `bench/bakeoff_chrom_ours.sh`
checks this itself (`git -C <repo> diff --quiet 3b096012 HEAD -- src Cargo.toml Cargo.lock`) before every
run and **refuses to run (exit 3)** if they differ. On every run it writes `run_provenance.txt` into the
label's own output directory (`$W/<label>/run_provenance.txt`) containing: the date; `git rev-parse HEAD` of
the repo; the src-identity check result; the binary's path and `sha256sum`; every currently-set `RUSTLE_*`
variable (after the unset step — so an unexpected non-empty listing here, beyond
`RUSTLE_LEGACY_PLACEMENT_DEDUP`, is itself evidence of a problem); and the exact command line used.

**Tool versions.** gffcompare: pinned by the environment to v0.12.10; `bench/chrom_score.sh` runs
`gffcompare --version` fresh for every label and writes it to `$W/gffcompare/<label>.provenance.txt`, along
with (when `--sqanti` is passed) `git -C <SQANTI3 checkout> rev-parse HEAD` of the SQANTI3 install actually
used — the specific commit is recorded there at run time rather than hardcoded here, since a conda-env
reinstall could change it between now and when chr17 is actually run; a differing recorded commit across
labels would itself be a reportable anomaly. gffcompare is invoked as `gffcompare -r <ref> -o <label> <gtf>`
with no other flags (no `-R`/`-Q` — either would change precision and confound E2). SQANTI3 is invoked with
exactly `--report skip -t 4` (plus `--isoforms`/`--refGTF`/`--refFasta`/`-o`/`-d`). StringTie and FLAIR run
exactly as in `bench/bakeoff_chr20_stringtie.sh` / `bench/bakeoff_chr20_flair.sh` (same flags, same FLAIR
`correct`-skip and PATH workarounds); their versions are not separately pinned or logged, since both arms are
context only and never enter a verdict.

**Verdict-script invocation and argument mapping** (copied from the plan's Task 12 Step 1,
`docs/superpowers/plans/2026-09-16-gtf-refine-and-dedup-fix.md`; gffcompare writes the `.tmap`/`.refmap` next
to the query GTF, i.e. in `$W/<label>/`, as `<label>.<gtf basename>.tmap`, NOT in `$W/gffcompare/`):
```
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17
BT=$(ls $W/baseline/baseline.*.tmap); NT=$(ls $W/bundle/bundle.*.tmap)
python3 bench/gtf_refine_verdict.py "$BT" $W/gffcompare/baseline.stats $W/sqanti3/baseline/baseline_junctions.txt \
  "$NT" $W/gffcompare/bundle.stats $W/sqanti3/bundle/bundle_junctions.txt \
  $W/sqanti3/baseline/baseline_classification.txt $W/sqanti3/abl_tss/abl_tss_classification.txt > /tmp/t9_verdict.json
```
Argument mapping: `BASE_TMAP` = baseline's `.tmap` (in `$W/baseline/`); `BASE_STATS` =
`$W/gffcompare/baseline.stats`; `BASE_JUNCTIONS` = `$W/sqanti3/baseline/baseline_junctions.txt`;
`BUNDLE_TMAP` = bundle's `.tmap` (in `$W/bundle/`); `BUNDLE_STATS` = `$W/gffcompare/bundle.stats`;
`BUNDLE_JUNCTIONS` = `$W/sqanti3/bundle/bundle_junctions.txt`; `BASE_CLASSIFICATION` =
`$W/sqanti3/baseline/baseline_classification.txt`; `ABL_TSS_CLASSIFICATION` =
`$W/sqanti3/abl_tss/abl_tss_classification.txt` (E5/TSS verdict only — BASELINE's classification file is
reused as the E5 baseline arm; ABL_TSS supplies the only-5'-end-differs comparison).

**Failure policy.** Infrastructure failures (a tool fails to install, a path is wrong, disk fills up, a
process is killed by the OOM/crash-avoidance rules) may be fixed and the run repeated — that is not a rule
change. In contrast: ANY change to `src/`, `Cargo.toml`, `Cargo.lock`, any of the five scripts above, or
`bench/gtf_refine_verdict.py` AFTER a chr17 output has been generated voids the test for every arm already
run — it is reported as VOID, not silently re-run and not quietly patched over. An E4 novel-junction
violation is scored as E4 FAILING regardless of its cause (implementation bug or otherwise) — a violation is
never explained away into a passing E4.

## Known limits (declared now)

- `tss` can only move a 5' end downstream (toward the read cluster) or leave it unchanged, never upstream of
  the most extreme exact read. On chr20, 22 of 348 chain-matched models had the reference TSS upstream of
  every exact read, so `tss` cannot reach the true TSS for them regardless of W (ceiling 326/348, not
  348/348) — `bench/CHR20_ASSEMBLER_COMPARISON.md`, "Known limits of `tss` (development observations, chr20,
  W = 50)" table (Task 7/9 section).
- `tss` can move a well-supported model AWAY from its true TSS when most reads are 5'-truncated: chr20
  example `DN_chr20_37420188_3` (132 reads, gene ROMO1) moved SQANTI3 `diff_to_TSS` from −8 to −168 bp — the
  majority of 5'-truncated reads outvoted the full-length minority. This is 1 of 6 `moved_out` cases at
  W = 50 on chr20 — same table.
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
