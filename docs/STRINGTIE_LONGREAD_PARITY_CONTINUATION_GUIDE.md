# StringTie Long-Read Parity Continuation Guide

This guide is a practical handoff for any AI agent (or human) to continue parity work between
Rustle and StringTie for **long-read de novo** and **guided** assembly.

It focuses on repeatable experiments, trace artifacts, and patch strategy rather than one-off fixes.

---

## 1) Objective

Improve Rustle so its long-read behavior converges toward StringTie (`stringtie -L`) by:

1. aligning intermediate decisions (junctions, cgroup partition, graph topology),
2. reducing unexplained divergences bundle-by-bundle,
3. validating with both instrumentation parity and final transcript comparisons.

Keep all changes measurable and reversible.

---

## 2) Current Parity Instrumentation (already available)

### StringTie (`stringtie/rlink.cpp`)

- `PARITY_PARTITION_TSV`
- `PARITY_JUNCTION_TSV`
- `PARITY_BUNDLE_ORDER_TSV`
- `PARITY_TRACE_TSV` + `PARITY_TRACE_LEVEL=1|2|3`
  - emits `trace_event_v1` rows (function/stage events)

### Rustle

- `RUSTLE_PARITY_PARTITION_TSV`
- `RUSTLE_PARITY_JUNCTION_TSV`
- `RUSTLE_PARITY_TRACE_TSV` + `RUSTLE_PARITY_TRACE_LEVEL=1|2|3`
  - emits:
    - `trace_event_v1`
    - `cgroup_summary_v1`
    - `graph_topology_v1` (`post_graph_build`, `pre_read_map`)
    - `read_map_summary_v1`
    - optional `read_exons_v1` (`RUSTLE_PARITY_READ_EXONS=1`)

### Analysis scripts (Rustle `scripts/`)

- `compare_partition_geometry.py`
- `parity_fidelity_report.py`
- `disambiguate_parity_bundles.py`

---

## 3) Important active parity knobs

These are key for current parity experiments:

- `RUSTLE_FREEZE_GRAPH_PRE_READ_MAP=1`
  - freezes coverage/terminal graph mutations before read mapping (Patch A behavior).
- `RUSTLE_KEEP_PRECANONICAL_JUNCTIONS=1`
  - uses pre-canonical junction stats as graph input source (Batch B2 probe).
- `RUSTLE_PARITY_FORCE_KEEP_JUNCTIONS_TSV=/path/file.tsv`
  - targeted keep list for specific bundle/junction tuples (Batch B2.2).
- `RUSTLE_PARITY_FORCE_DROP_JUNCTIONS_TSV=/path/file.tsv`
  - targeted drop list for specific bundle/junction tuples (Batch B2.2).
- `RUSTLE_CGROUP_STRICT_BOUNDARIES=1`
  - C4 probe: disables pair-stitch boundary extension in cgroup group construction.
  - currently keep this OFF by default (observed no partition gain and mild junction regression).
- `RUSTLE_STRINGTIE_PAIR_STITCH_COORDS=1`
  - C5 probe: use StringTie-style first-segment coordinates during pair-stitch boundary updates.
  - currently keep this OFF by default (no measured parity improvement).
- `RUSTLE_STRINGTIE_LONGREAD_KEEP_EXON=1`
  - C6 probe: enable StringTie long-read small-exon keep/drop clause in merge-to-existing-group path.
  - currently keep this OFF by default (measured regression on partition exact).
- `RUSTLE_STRINGTIE_MERGE_CLOSE_PASS=1`
  - C7 probe: use StringTie-style single-pass adjacency merge order in `merge_close_groups`.
  - currently keep this OFF by default (no measured parity improvement yet).
- `RUSTLE_MERGE_HARD_BOUNDARY_VETO=1`
  - C8 probe: block merge-close joins when left group has `hardend` or right group has `hardstart`.
  - currently keep this OFF by default (no measured parity improvement yet).
- `RUSTLE_GROUP_START_EXTEND_FIRST_EXON_ONLY=1`
  - C9 probe: only allow extending an existing cgroup start from exon 0 of the read.
  - currently keep this OFF by default (measured partition regression).

For keep/drop TSVs, columns are:
`chrom<TAB>bundle_start<TAB>bundle_end<TAB>donor<TAB>acceptor`
using Rustle internal junction coordinates.

Use these for controlled A/B runs; do not assume they are production defaults.

---

## 4) Canonical full-run recipe (GGO_19)

Run from repository root paths shown here.

```bash
# Paths
ST=/mnt/c/Users/jfris/Desktop/stringtie/stringtie
RU=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
BAM=/mnt/c/Users/jfris/Desktop/GGO_19.bam
OUT=/mnt/c/Users/jfris/Desktop/ggo19_trace_full_run
mkdir -p "$OUT"

# 1) StringTie
PARITY_TRACE_TSV="$OUT/st_trace.tsv" \
PARITY_TRACE_LEVEL=2 \
PARITY_PARTITION_TSV="$OUT/st_part.tsv" \
PARITY_JUNCTION_TSV="$OUT/st_junc.tsv" \
"$ST" -L "$BAM" -o "$OUT/st.gtf"

# 2) Rustle (baseline parity run with current best probes)
RUSTLE_PARITY_TRACE_TSV="$OUT/ru_trace.tsv" \
RUSTLE_PARITY_TRACE_LEVEL=2 \
RUSTLE_PARITY_PARTITION_TSV="$OUT/ru_part.tsv" \
RUSTLE_PARITY_JUNCTION_TSV="$OUT/ru_junc.tsv" \
RUSTLE_FREEZE_GRAPH_PRE_READ_MAP=1 \
RUSTLE_KEEP_PRECANONICAL_JUNCTIONS=1 \
"$RU" -L "$BAM" -o "$OUT/ru.gtf"

# 3) Summary report
python3 /mnt/c/Users/jfris/Desktop/Rustle/scripts/parity_fidelity_report.py \
  --partition-st "$OUT/st_part.tsv" --partition-ru "$OUT/ru_part.tsv" \
  --junction-st "$OUT/st_junc.tsv" --junction-ru "$OUT/ru_junc.tsv" \
  > "$OUT/parity_report.txt"

# 4) Bundle-level disambiguation
python3 /mnt/c/Users/jfris/Desktop/Rustle/scripts/disambiguate_parity_bundles.py \
  --partition-st "$OUT/st_part.tsv" --partition-ru "$OUT/ru_part.tsv" \
  --junction-st "$OUT/st_junc.tsv" --junction-ru "$OUT/ru_junc.tsv" \
  --stringtie-trace "$OUT/st_trace.tsv" --rustle-trace "$OUT/ru_trace.tsv" \
  -o "$OUT/disambig.tsv"
```

### Focused B2.2 run template (targeted keep/drop)

```bash
RUSTLE_FREEZE_GRAPH_PRE_READ_MAP=1 \
RUSTLE_KEEP_PRECANONICAL_JUNCTIONS=1 \
RUSTLE_PARITY_FORCE_KEEP_JUNCTIONS_TSV="$OUT/force_keep.tsv" \
RUSTLE_PARITY_FORCE_DROP_JUNCTIONS_TSV="$OUT/force_drop.tsv" \
RUSTLE_PARITY_TRACE_TSV="$OUT/ru_trace_batch.tsv" \
RUSTLE_PARITY_TRACE_LEVEL=2 \
RUSTLE_PARITY_PARTITION_TSV="$OUT/ru_part_batch.tsv" \
RUSTLE_PARITY_JUNCTION_TSV="$OUT/ru_junc_batch.tsv" \
"$RU" -L "$BAM" -o "$OUT/ru_batch.gtf"
```

---

## 5) How to continue (recommended loop)

Use this 6-step loop for each batch:

1. **Rank offenders**
   - Use `disambig.tsv`.
   - Prioritize loci with:
     - `junction_slop2_bipartite_ok=0` (true junction divergence),
     - high `ru_read_mapped_tfs_sum`,
     - graph change flags (`ru_graph_topology_changed_pre_read_map=1`).

2. **Localize drop stage**
   - For each missing ST junction, classify first stage where Rustle drops it:
     - before `read_union`,
     - `read_union -> pre_canonical`,
     - `pre_canonical -> graph_input`,
     - `graph_input -> cgroup_feed`.

3. **Patch smallest shared mechanism**
   - Prefer mechanism-level fixes over per-locus hacks.
   - Keep patch scope narrow and env-gated when uncertain.

4. **Run full GGO_19 parity**
   - Re-generate all TSVs and disambiguation.
   - Confirm targeted metric improvement and watch for regressions.

5. **Check transcript-level impact**
   - Compare Rustle vs StringTie GTF using `gffcompare`.
   - Ensure parity changes do not collapse transcript quality.

6. **Promote, gate, or revert**
   - If improvement is local and safe: keep.
   - If mixed: keep behind env flag.
   - If broad regression: revert quickly.

---

## 6) Current known residual pattern

After Patch A + B1 + B2 + focused B2.2 keep/drop probes, residual hard divergences are small and localized.

### Latest measured snapshot (full `GGO_19.bam`)

- `partition exact`: **784 / 1000 (78.4%)**
- best stage-pair junction exact (slop0): **97.7%**
  - `ST aggregate_graph_ready_after_read_pass` vs `RU cgroup_feed`
- mean Jaccard at that pair: **99.5%**
- hard bundle mismatches (`junction_slop2_bipartite_ok=0`): **23**

Compared with earlier checkpoints in this session:

- pre-B2.2 focused work: hard mismatches ~33
- after Patch A + B1 + B2: hard mismatches ~28
- after B2.2 targeted keep/drop: hard mismatches ~23

Typical remaining causes:

- junctions dropped at `read_union -> pre_canonical`,
- junctions dropped between `graph_input` and `cgroup_feed`,
- occasional read-correction-level miss before `read_union`.

Three practical unresolved classes remain:

1. ST-only junction missing in Rustle cgroup feed at a few loci.
2. Rustle-only extra junction retained in cgroup feed at a few loci.
3. Chain/topology partition differences where junction tables already match.

Interpretation: this is no longer broad coordinate-convention noise; now it is specific stage logic.

---

## 7) De novo vs guided parity strategy

### De novo (`-L` without guide)

Primary target first. Validate:

- partition/junction parity,
- event-stage consistency,
- `gffcompare` transcript sensitivity/precision vs StringTie output/reference.

### Guided (`-L -G`)

After de novo stabilizes:

1. add paired guided runs with same parity dumps,
2. include guide-related checkpoints in trace events (guide acceptance, guide boundary split usage),
3. compare guide-sensitive divergences separately from de novo core logic.

Do not mix de novo and guided conclusions without labeling.

---

## 8) Suggested next experiments (high value)

1. **Residual hard-junction micro-panel**
   - Build mini-BAMs for top unresolved loci.
   - Iterate faster than full BAM for mechanism testing.

2. **Stage-specific forced keep/drop probes**
   - Use temporary env-gated keep/drop maps for diagnostics only.
   - Replace with generalized logic once pattern repeats.

3. **StringTie trace enrichment**
   - Add event rows around any stage where Rustle still cannot map causality one-to-one.

4. **Guided parity extension**
   - Add and compare guide-event rows in both tools.

### Batch C kickoff note (post-B2.2d)

- Generated a fresh C queue:
  - `ggo19_trace_full_run/batchC_top23.tsv`
  - `ggo19_trace_full_run/batchC_top23.md`
- First C probe added:
  - `RUSTLE_PARTITION_SIGNATURE_MODE=subbundle|color_components`
  - location: parity partition signature selection in `pipeline.rs`
- Initial result on full `GGO_19` with current B2.2d knobs:
  - switching to `color_components` did **not** change headline partition exact
    (stayed `784/1000`, 78.4%).
- Interpretation: remaining C gaps are likely not a simple signature-encoding mode toggle;
  continue with cgroup/subbundle construction logic and per-bundle chain segmentation.

### Batch C2 result (chain-detail trace)

- Added `partition_chain_detail_v1` rows in Rustle trace (`RUSTLE_PARITY_CHAIN_DETAIL=1`):
  - mode `subbundle`
  - mode `color_components`
- Generated:
  - `ggo19_trace_full_run/batchC_top23.tsv`
  - `ggo19_trace_full_run/batchC2_chain_report.md`
- Finding on top-23 C loci:
  - ST chain count == RU subbundle chain count == RU color-component chain count in all 23.
  - Therefore remaining C mismatches are likely **segment boundary/grouping details inside chains**,
    not chain-count/cardinality.
- Metrics stayed stable during C2 instrumentation run:
  - hard mismatches (`junction_slop2=0`): 23
  - junction exact count: 193/216 mismatch bundles.

### Batch C3 result (segment drift quantification)

- Added script:
  - `scripts/analyze_partition_segment_drift.py`
- Outputs from latest run:
  - `ggo19_trace_full_run/batchC3_segment_drift.tsv`
  - `ggo19_trace_full_run/batchC3_segment_drift.md`
- On Batch C top-23 (junction-exact, partition-mismatch loci):
  - all 23 have equal ST/RU segment counts per bundle.
  - dominant RU-ST start delta: `0` (vast majority), with small negative tails (`-1..-5`).
  - dominant RU-ST end delta: `0` (vast majority), with small positive tails (`+1..+5`).
- Interpretation:
  - Batch C residual is mainly **small boundary expansion/contraction inside matched chain counts**.
  - Next C patch should target boundary placement rules in cgroup/bundlenode construction, not
    chain cardinality or signature mode selection.

### Batch C4 result (strict-boundary probe)

- Added env-gated probe in `bundle_builder.rs`:
  - `RUSTLE_CGROUP_STRICT_BOUNDARIES=1`
  - disables two pair-fragment stitching boundary-extension paths that can widen group ends.
- Full `GGO_19` run artifacts:
  - `ggo19_trace_full_run/ru_part_batchC4.tsv`
  - `ggo19_trace_full_run/ru_junc_batchC4.tsv`
  - `ggo19_trace_full_run/ru_trace_batchC4.tsv`
  - `ggo19_trace_full_run/parity_report_batchC4.txt`
  - `ggo19_trace_full_run/disambig_batchC4.tsv`
  - `ggo19_trace_full_run/batchC4_segment_drift.tsv`
- Result:
  - partition exact stayed flat at `784/1000` (`78.4%`).
  - best junction exact at the best stage pair remained `97.2%` (no improvement).
  - hard mismatches (`junction_slop2=0`) increased to `28` (from prior `23`).
  - junction-exact+partition-mismatch queue moved to `188` bundles (from `193`) with no net headline gain.
- Decision:
  - keep C4 probe env-gated and OFF for baseline parity runs.
  - next C iteration should prefer **local/conditional boundary normalization** keyed by evidence
    (e.g., cgroup-specific support) instead of global pair-stitch suppression.

### Batch C5 result (StringTie pair-stitch coordinate probe)

- Added env-gated probe in `bundle_builder.rs`:
  - `RUSTLE_STRINGTIE_PAIR_STITCH_COORDS=1`
  - uses first-exon coordinates for pair-stitch extension points, to mirror StringTie semantics.
- Full `GGO_19` run artifacts:
  - `ggo19_trace_full_run/ru_part_batchC5.tsv`
  - `ggo19_trace_full_run/ru_junc_batchC5.tsv`
  - `ggo19_trace_full_run/ru_trace_batchC5.tsv`
  - `ggo19_trace_full_run/parity_report_batchC5.txt`
  - `ggo19_trace_full_run/disambig_batchC5.tsv`
- Result:
  - partition exact stayed flat at `784/1000` (`78.4%`).
  - best junction exact remained `97.2%`.
  - hard mismatches (`junction_slop2=0`) remained elevated at `28`.
- Decision:
  - keep C5 env-gated and OFF for baseline parity runs.
  - prioritize next C patches in `merge_close_groups` and cgroup overlap/keep-exon guards,
    where many same-count topology mismatches still originate.

### Batch C6 result (StringTie long-read keep-exon clause probe)

- Added env-gated probe in `bundle_builder.rs`:
  - `RUSTLE_STRINGTIE_LONGREAD_KEEP_EXON=1`
  - enables long-read `CHI_THR`/`DROP` small-exon filter in the merge-to-existing-group keep-exon path.
- Full `GGO_19` run artifacts:
  - `ggo19_trace_full_run/ru_part_batchC6.tsv`
  - `ggo19_trace_full_run/ru_junc_batchC6.tsv`
  - `ggo19_trace_full_run/ru_trace_batchC6.tsv`
  - `ggo19_trace_full_run/parity_report_batchC6.txt`
  - `ggo19_trace_full_run/disambig_batchC6.tsv`
- Result:
  - partition exact regressed to `762/1000` (`76.2%`) from `78.4%`.
  - best junction exact stayed `97.2%` (no junction gain).
- Decision:
  - keep C6 env-gated and OFF for baseline parity runs.
  - this suggests Rustle’s current long-read keep behavior in that path was closer to observed parity than
    a direct activation of this StringTie clause; next patch should focus on merge-boundary semantics instead.

### Batch C7 result (StringTie merge-close pass order)

- Added env-gated probe in `bundle_builder.rs`:
  - `RUSTLE_STRINGTIE_MERGE_CLOSE_PASS=1`
  - switches merge pass to StringTie-style single-pass adjacent `last/proc` iteration order.
- Full `GGO_19` run artifacts:
  - `ggo19_trace_full_run/ru_part_batchC7.tsv`
  - `ggo19_trace_full_run/ru_junc_batchC7.tsv`
  - `ggo19_trace_full_run/ru_trace_batchC7.tsv`
  - `ggo19_trace_full_run/parity_report_batchC7.txt`
  - `ggo19_trace_full_run/disambig_batchC7.tsv`
- Result:
  - partition exact stayed flat at `784/1000` (`78.4%`).
  - best junction exact stayed `97.2%`.
  - hard mismatches remained at `28`.
- Decision:
  - keep C7 env-gated and OFF for baseline parity runs.
  - next candidate should target missing guide-boundary merge blockers (boundary_left/right semantics)
    or more local boundary writes in cgroup construction.

### Batch C8 result (hard-boundary merge veto)

- Added env-gated probe in `bundle_builder.rs`:
  - `RUSTLE_MERGE_HARD_BOUNDARY_VETO=1`
  - prevents close-group merge when merge crosses an explicit hard boundary marker.
- Full `GGO_19` run artifacts:
  - `ggo19_trace_full_run/ru_part_batchC8.tsv`
  - `ggo19_trace_full_run/ru_junc_batchC8.tsv`
  - `ggo19_trace_full_run/ru_trace_batchC8.tsv`
  - `ggo19_trace_full_run/parity_report_batchC8.txt`
  - `ggo19_trace_full_run/disambig_batchC8.tsv`
- Result:
  - partition exact stayed flat at `784/1000` (`78.4%`).
  - best junction exact stayed `97.2%`.
  - hard mismatches stayed `28`.
- Decision:
  - keep C8 env-gated and OFF for baseline parity runs.
  - this indicates merge-close boundary veto alone is insufficient; next step should target
    per-read group-boundary updates in `process_read_for_group` with event-level first-divergence tracing.

### Batch C9 result (first-exon-only start extension)

- Added env-gated probe in `bundle_builder.rs`:
  - `RUSTLE_GROUP_START_EXTEND_FIRST_EXON_ONLY=1`
  - restricts `group.start` extension in overlap path to first exon only.
- Full `GGO_19` run artifacts:
  - `ggo19_trace_full_run/ru_part_batchC9.tsv`
  - `ggo19_trace_full_run/ru_junc_batchC9.tsv`
  - `ggo19_trace_full_run/ru_trace_batchC9.tsv`
  - `ggo19_trace_full_run/parity_report_batchC9.txt`
  - `ggo19_trace_full_run/disambig_batchC9.tsv`
- Result:
  - partition exact regressed to `777/1000` (`77.7%`) from `78.4%`.
  - best junction exact remained `97.2%`.
  - hard mismatches remained `28`.
- Decision:
  - keep C9 env-gated and OFF for baseline parity runs.
  - this suggests non-first-exon start updates are needed for current parity; next probes should
    focus on selective conditions (junction-support / pair-status) rather than a global first-exon restriction.

---

## Batch D: Pure-source-stub gate experiments (2026-04-23)

**Baseline state at start:** partition_exact=78.4%, junction_exact=97.2%, Rustle 1914 transcripts,
StringTie 1839 transcripts, gffcompare Sn/Pr=89.0/85.5 (transcript level).

### False-positive analysis

gffcompare class codes for non-matching Rustle transcripts (277 total):
- `j` partial intron chain: 175 (flow=114, checktrf=46, orphan=15)
- `k` contained w/ diff junction count: 30 (flow=27, checktrf=3)
- `m` retained intron: 21 (flow=13, checktrf=8)
- `x` opposite strand: 19 (ALL flow) ← cleanest target; includes RSTL.342.1 and RSTL.40 cluster
- `c` contained: 14 (checktrf=8, orphan=4, flow=2)
- `n` intronic: 12 (flow=7, checktrf=5)
- `u`/`o` other: 6

The 19 `x`-class (wrong-strand) flow transcripts are diagnostically interesting:
- RSTL.40 cluster: 5 plus-strand transcripts at 20567313-20581018 (StringTie has ONLY minus-strand here)
- RSTL.342.1: 1 minus-strand at 22153090-22155217 (StringTie has plus-strand)
- Pattern: reads in regions between two same-strand genes get assigned to the wrong strand in Rustle's
  sub-bundle context; in StringTie's full-bundle context, they'd be absorbed into the dominant strand.

### Experiment D1: `RUSTLE_BLOCK_ARTIFICIAL_SOURCE_STUB=1` (back stub → checktrf)

- Condition: when `back_used_pure_source_stub=true` AND first real node has no hardstart/hardend,
  block direct-store and route to checktrf.
- Result: Rustle 1952 transcripts (flow=837, checktrf_rescue=1093 ← UP from 260!)
- Problem: checktrf RESCUES the blocked seeds instead of matching them to existing flow transcripts,
  because in isolated sub-bundles there are no matching transcripts to absorb them.
- Decision: **REVERTED**

### Experiment D2: `RUSTLE_DROP_BOTH_ARTIFICIAL_STUBS=1` (both stubs → drop, no path limit)

- Condition: when BOTH `back_used_pure_source_stub=true` AND `fwd_used_pure_sink_stub=true` AND
  no hardstart/hardend, DROP seed entirely (skip checktrf too).
- Result: Rustle 1213 transcripts; gffcompare Sn=49.1, Pr=74.4 → MUCH WORSE
- Problem: `fwd_used_pure_sink_stub=true` fires for almost ALL seeds (any seed whose last fwd step
  uses a 2-node [last_real_node → sink] stub), so the "both stubs" condition is very broad.
- Decision: **REVERTED**

### Experiment D3: `RUSTLE_DROP_BOTH_ARTIFICIAL_STUBS=1` (both stubs AND path.len() ≤ 4)

- Added `path.len() <= 4` guard (targets only 1-2 real node transcripts spanning whole sub-bundle).
- Result: Rustle 1829 transcripts (close to StringTie 1839);
  gffcompare Sn=84.6, Pr=85.1 → Sn regressed -4.4 points (1637 → 1556 matching chains)
- Problem: path.len() ≤ 4 still drops many legitimate 2-exon transcripts that happen to span
  their entire sub-bundle. Cannot distinguish artificial truncation from genuine short transcripts.
- Decision: **REVERTED**

### Root cause conclusion

The path_extract gate approach (at the back/fwd certification step) cannot fix sub-bundle inflation
because:
1. **back/fwd stub conditions fire broadly** — most seeds in sub-bundles near boundaries use stubs
2. **checktrf rescues instead of matches** — isolated sub-bundles have no matching transcripts, so
   checktrf creates new independent transcripts instead of absorbing blocked seeds
3. **stub conditions don't distinguish artificial vs. legitimate** — path.len() is not a reliable
   discriminator; many legitimate 2-exon transcripts span their sub-bundle exactly

**The true root cause is architectural: sub-bundle splitting creates isolated processing contexts
with artificial source/sink nodes.** Seeds in these contexts trivially certify and get stored.
This produces both:
- Truncated transcripts (reads with exons outside the sub-bundle)
- Wrong-strand transcripts (reads in sparse regions between same-strand genes)

### Experiment D4: `RUSTLE_NO_BAD_JUNCTION_COLOR_SPLIT=1` (suppress bad-junc color splits)

Root cause investigation: sub-bundle count is driven by "bad junction" color splitting in:
- `bundle.rs` (non-3strand path): lines 764 (adopt without merge) and 795 (new color on no-overlap)
- `bundle_builder.rs` `process_read` function: line ~1496 (new color on bad-junction no-overlap)

StringTie DOES NOT split on bad junctions during bundle formation. It uses overlap-based merging.

Gate: suppress bad-junction fresh color creation in both files.

Result:
- Partition exact: 69.2% (DOWN from 78.4%, -9.2 points)
- gffcompare: Sn=89.0/Pr=85.1 (nearly identical to baseline 89.0/85.5)
- Matching transcripts: 1637 (SAME as baseline)
- Transcript count: 1923 (+9 from baseline 1914)

**Why it regressed partition:** Merging sub-bundles creates larger units whose internal bundlenode
structure doesn't match StringTie's partition. StringTie's 2 cgroups for 743 kb loci have specific
internal geometry from their own algorithm; Rustle's merged 2-sub-bundle has different geometry.
Even though sub-bundle count decreased, the local match quality dropped.

**Key insight:** The bad-junction color split IS actually helping partition geometry by creating
small, locally-correct sub-bundles that more closely match StringTie's local structure. Removing
these splits creates larger sub-bundles with different geometry.

Decision: **REVERTED**

### Recommended next approach

All path_extract gates and bundle-formation gates tested so far have failed or regressed.

Deeper root cause:
- Top offenders: StringTie uses 2 cgroups for loci where Rustle creates 60-90 sub-bundles
  (disambig `ru_cgroup_sub_bundles` column; e.g., locus 40036006-40779205 has 90 sub-bundles)
- The partition exact measures LOCAL bundlenode geometry, which depends on fine-grained sub-bundle
  structure. Merging sub-bundles changes this local structure.
- StringTie's 2 cgroups correspond to its internal algorithm; Rustle can't simply replicate this
  by merging sub-bundles.

Consider:
1. **Understand StringTie's cgroup algorithm more deeply**: What makes StringTie create exactly 2
   cgroups for a 743kb locus? How does it process reads differently from Rustle?
   StringTie likely uses strand-based grouping with aggressive overlap-based merging.

2. **Cross-boundary read truncation detection**: When a read's CIGAR extends beyond the sub-bundle
   boundaries, flag its derived transfrag. Do not direct-store transcripts derived entirely from
   cross-boundary reads (their exon pattern is artificially truncated). This could target the
   specific `x`-class and truncated `j`-class false positives without affecting partition geometry.

3. **Wrong-strand isolation fix**: When a sub-bundle has low abundance AND is between two same-strand
   bundles of the opposite strand, absorb its reads into the dominant-strand context rather than
   creating an independent opposite-strand sub-bundle.

4. **Accept the current parity baseline**: Sn/Pr=89.0/85.5 with partition 78.4% is already good.
   The remaining 216 bundle mismatches are rooted in the fundamental architectural difference between
   StringTie's single-bundle processing and Rustle's sub-bundle approach. Without deeper changes to
   the core path extraction algorithm (not just the bundling), further gains may be incremental.

---

## 9) Guardrails for any agent

- Always run `cargo build --release` after Rustle code changes.
- Prefer env-gated experimental behavior until metrics are stable.
- Keep one change per batch so attribution is clear.
- Save outputs to unique files (`ru_*_batchX.tsv`) to preserve comparability.
- Never claim parity progress without both:
  - instrumentation metrics (`parity_fidelity_report.py`),
  - transcript-level check (`gffcompare`).

---

## 10) Quick status checklist template

For each batch, record:

- Change:
- Env flags used:
- `partition exact %`:
- `best junction exact % (slop0)`:
- `hard mismatches (junction_slop2=0)`:
- top 5 offender changes:
- `gffcompare` transcript-level change:
- keep/revert decision:

Use this file as the canonical handoff update log.

