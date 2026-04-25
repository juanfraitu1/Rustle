# Sub-Bundle Architectural Parity — Handoff (sessions 1–7, 2026-04-24)

This document is the practical handoff for any AI agent (or human) continuing
the long-read parity investigation between Rustle and StringTie (`-L` mode).
It captures the architectural patterns that DID and DIDN'T work across seven
investigation sessions.

**Read first**: `docs/STRINGTIE_LONGREAD_PARITY_CONTINUATION_GUIDE.md` for
broader context. This document is the focused architectural addendum.

---

## 1. Current state

On `GGO_19.bam` vs `stringtie_parity_v2.gtf`:

| metric | start of investigation | end of session 7 | Δ |
|---|---|---|---|
| Matching intron chains | 1633 | **1662** | +29 |
| Intron chain Sn | 88.8 | **91.6** | +2.8 |
| Intron chain Pr | 85.4 | **86.0** | +0.6 |

To run with all session 7 wins:

```bash
RUSTLE_ALT_SPLICE_RESCUE=1 \
  /path/to/rustle -L -o out.gtf input.bam
```

To restore the pre-session-7 baseline (1655/91.2/85.2):

```bash
RUSTLE_LEGACY_SUB_JUNC_STATS=1 \
RUSTLE_LEGACY_SUB_REDIRECT_MAP=1 \
RUSTLE_SUB_GOOD_JUNC_SECOND_STAGE=1 \
RUSTLE_SUB_JUNC_DEMOTIONS=1 \
RUSTLE_SUB_REPAIR_READS=1 \
RUSTLE_SUB_CORRECT_JUNCTIONS=1 \
RUSTLE_SUB_SNAP_GUIDES=1 \
RUSTLE_SUB_APPLY_FILTERS=1 \
RUSTLE_ALT_SPLICE_RESCUE=1 \
  /path/to/rustle -L -o out.gtf input.bam
```

Each opt-out flag is independent — toggle one at a time to A/B test individual
architectural changes.

---

## 2. The single insight that unlocked everything

After **eight regressing port attempts** (shared-site-weak filter,
good-junc-strict, longtrim-strict-lend, badjunc-trust-node-boundary,
dual-junction-set, alt-acceptor-rescue, BUNDLE_GRAPH disable, EXON_SKIP_OFF),
session 7 found the architectural root cause:

> **Rustle re-filters junction stats per sub-bundle from a narrower read pool.
> StringTie processes everything at the bundle level.**

When sub-bundles recompute junction stats:
- A junction supported by 100 reads in the outer bundle but only 1 read in
  this specific sub-bundle gets `nreads_good < threshold` → marked killed.
- The kill cascades: sub-bundle graph misses this junction → no node
  boundary at the donor → reads using this junction map to a different node
  path → transfrag dedup fragments → parse_trflong can't reconstruct the ST
  transcript.

This affected **74% of reads** in some bundles (e.g., the STRG.18 case at
`NC_073243.2:20126937-20733034 +`).

### The fix pattern

**Use outer-bundle junction stats throughout sub-bundle processing**, with
narrow exceptions where sub-bundle context is genuinely needed.

For each sub-bundle pass:
1. Identify whether outer bundle already does this work
2. If yes, gate the sub-bundle re-run as a no-op
3. Validate: `RUSTLE_LEGACY_<NAME>=1` should restore the pre-change behavior
   exactly

This worked for 8 sub-bundle passes (see §4 below). Only `count_good_junctions`
genuinely needs sub-bundle context (DEL-aware per-read resolution that
mutates the per-read state inside the sub-bundle).

---

## 3. Patterns that DON'T work (8 negative results)

When you find a sub-bundle decision that diverges from outer, do **NOT**
implement these:

| approach | why it fails |
|---|---|
| Inject killed junctions back at filter point | Cascading flow-graph effects |
| Lower threshold gates (min_reads, ratio, window) | Picks up too much noise |
| Disable Rustle's compensating filters | Rustle's tuning depends on them |
| Trust node boundary over kill flag | Reads then traverse non-existent edges |
| Strong-kills-weak shared-splice-site filter | Kills legitimate alt sites |

The pattern: **rescuing or relaxing AFTER sub-bundle kills happened produces
cascading regressions**. The kills themselves are the problem.

The thing that actually works: **don't kill in sub-bundle in the first
place**. Use outer-bundle decisions throughout.

---

## 4. The 8 default-on architectural alignments

Each is a sub-bundle pass that's now skipped or re-seeded from outer.

| # | Sub-bundle pass | Default behavior | Opt-out flag | Δ matches |
|---|---|---|---|---|
| 1 | `compute_initial_junction_stats_for_reads` | Seed from outer restricted | `RUSTLE_LEGACY_SUB_JUNC_STATS=1` | +2 |
| 2 | Final `sub_junction_stats` filter | Replace with outer restricted | (same flag) | (paired with #1) |
| 3 | `apply_junction_filters_and_canonicalize` | Skip (output discarded) | `RUSTLE_SUB_APPLY_FILTERS=1` | 0 |
| 4 | `correct_junctions_with_map` | Skip | `RUSTLE_SUB_CORRECT_JUNCTIONS=1` | 0 |
| 5 | `snap_junctions_to_guides` | Skip | `RUSTLE_SUB_SNAP_GUIDES=1` | 0 |
| 6 | `good_junc_second_stage` | Skip | `RUSTLE_SUB_GOOD_JUNC_SECOND_STAGE=1` | +1 |
| 7 | `apply_higherr_demotions`+`demote_runthrough`+`good_junc` | Skip | `RUSTLE_SUB_JUNC_DEMOTIONS=1` | +2 |
| 8 | `repair_reads_after_junction_quality` | Skip | `RUSTLE_SUB_REPAIR_READS=1` | 0 (Pr +0.3) |
| 9 | `junction_redirect_map` | Inherit outer | `RUSTLE_LEGACY_SUB_REDIRECT_MAP=1` | 0 |

Cumulative: +5 chains, +0.3 Sn, +0.6 Pr (all from `compute_initial`+`final
filter`, the rest are no-op simplifications or precision improvements).

`count_good_junctions` is intentionally NOT in this list. Skipping it
regresses -2 matches because it does DEL-aware per-read junction resolution
on the sub-bundle's specific read context.

---

## 5. Where to look next (ranked by promise)

### 5.1 Long_max_flow numerical parity

ST's `long_max_flow` and Rustle's `long_max_flow_direct` are mathematically
equivalent (same Edmonds-Karp augmenting path, same coverage formula). But
seed-order-dependent nodecov depletion produces different per-seed flux
values. See `project_rustle_parity.md` "Investigation into flow coverage"
section.

**Concrete starting point**: bundle 16472610-16819296 minus, STRG.18.1
case (cov 4.52 in ST, ~2.50 in Rustle after depletion). Existing dumps:

- ST: `PARITY_NODEFLUX_DUMP=/path/file.tsv stringtie -L ...`
- Rustle: `RUSTLE_PATH_DECISION_TSV=/path/file.tsv rustle -L ...`

Diff seed-by-seed using coordinate-based identity (chrom, ref_start,
ref_end, n_nodes). A specific divergence is documented at the STRG.18 case.

**Why I didn't pursue**: the per-seed depletion ordering is sensitive to
parse_trflong's seed iteration order, which depends on hash-map iteration
in process_transfrags. Aligning is a multi-session port of a substantial
helper.

### 5.2 BAM read ingest differences

ST and Rustle differ slightly in which reads they keep:

- ST: 5228 reads in bundle 16472610-16819296
- Rustle: ~5200 (varies by version)

Differences come from:
- Secondary/supplementary alignment handling
- Pair-aware processing
- Strand inference fallback (polyA/polyT tail evidence)

ST's instrumentation: `PARITY_BUNDLE_TRACE=/path/file.tsv` traces every
read processed in the reader loop with bundle-state info.

**Concrete approach**: dump the read set for each bundle from both tools,
diff by qname (or coord+CIGAR signature), look for systematic patterns.

### 5.3 Read→transfrag path divergence (the 4587 cases)

Even with all session-7 changes, ~4587 reads have different node paths
than ST. Documented breakdown in memory file:

- 5096 last_node_end drift (Rustle's last node extends past ST's by 20-100bp)
- 4587 rustle_splits_more (Rustle fragments reads more than ST)
- 1776 node_count_diff (ST has extra internal nodes)

Each has documented attempts and results in the memory file.

**Lower-risk subset**: `node_count_diff` cases where ST has 4 sub-nodes
where Rustle has 1 (the STRG.18 wobble pattern). With session-7 outer-stats
inheritance, some of these may now resolve naturally. Re-run the diff to
quantify how many remain.

### 5.4 Guide-driven runs (-G)

Session 7 only validated on de-novo (`-L` without `-G`). The
`mark_guide_junctions_for_junctions` and snap-to-guides paths were
no-ops in our test. For guide-driven runs:

- Verify outer's `mark_guide_junctions_for_junctions` flows correctly
  into sub-bundle stats via the new outer-stats inheritance.
- Verify guide-snap redirects in outer's `junction_redirect_map`
  propagate to sub-bundle via `RUSTLE_LEGACY_SUB_REDIRECT_MAP=0` default.

If guided runs regress, gate the relevant sub-bundle pass behind a
specific flag (e.g., `RUSTLE_SUB_MARK_GUIDES=1` for guided mode).

---

## 6. Investigation methodology that works

### 6.1 Before changing anything

1. Reproduce the current baseline metric (`gffcompare` matching chains).
2. Identify a specific bundle where parity diverges. Use existing dumps:
   - `PARITY_PATH_DECISION_TSV` / `RUSTLE_PATH_DECISION_TSV` (per-seed)
   - `PARITY_GRAPH_DUMP` / `RUSTLE_PARITY_GRAPH_TSV` (graph nodes)
   - `PARITY_TF_DUMP` / `RUSTLE_PARITY_TF_TSV` (transfrags)
   - `PARITY_READ_NODE_TSV` / `RUSTLE_READ_NODE_TSV` (read→node paths)
   - `PARITY_JUNCTION_TSV` / `RUSTLE_PARITY_JUNCTION_TSV` (per-stage junctions)
   - `PARITY_CBUNDLE_DUMP` (ST CBundle structure)
3. Diff coordinate-normalize: ST 1-based inclusive ↔ Rustle 0-based half-open
   (start: ST_start - 1 = Rustle_start; end: same).

### 6.2 When investigating a divergence

Walk the bundle through each phase. Phases (from ST `build_graphs`):
1. Bundle setup
2. CPAS clustering (long reads)
3. Guide marking
4. Junction sweep (leftsupport/rightsupport aggregation, dinucleotide consensus)
5. Junction filter (`good_junc`)
6. Per-read changeleft/changeright redirect
7. `add_read_to_group` (color union-find)
8. `merge_fwd_groups`
9. Cross-strand color projection → CBundles
10. `create_graph` per CBundle (with longtrim inline)
11. `get_read_to_transfrag` (per-read → node-path)
12. `process_refguides` → `process_transfrags`
13. `parse_trflong` (seed iteration, prediction)

At each phase, dump state from both tools, diff, find the earliest
divergence. Fix THAT phase, not downstream.

### 6.3 When implementing a fix

1. **Always gate behind an env flag.** Use `RUSTLE_<DESCRIPTIVE_NAME>=1`
   as the opt-IN, or `RUSTLE_LEGACY_<NAME>=1` as the opt-OUT (when
   default-on).
2. **Verify the legacy flag exactly restores prior behavior.** Run twice:
   once with default flag, once with legacy flag. Metrics must match
   prior baseline byte-for-byte.
3. **Don't promote to default-on until verified across runs.**
4. **Document the flag in this guide** with the matches Δ and rationale.

### 6.4 Red flags — back out immediately

- Sn drops by more than +0 chains (any direct port that loses sensitivity
  is wrong even if Pr improves; precision is easy, sensitivity is hard).
- Test with all flags on legacy: must give identical pre-change baseline.
- Per-bundle diff worsens (use the read→node-path diff script).

---

## 7. Tools and scripts

### Building

```bash
cd /path/to/Rustle && cargo build --release
```

Build takes ~90s. Always run from inside `Rustle/` directory.

### Quick metric check

```bash
RUSTLE_ALT_SPLICE_RESCUE=1 \
  /path/to/rustle/target/release/rustle -L -o /tmp/ru.gtf /path/to/in.bam

cd /tmp && gffcompare -r /path/to/stringtie_parity.gtf -o cmp ru.gtf
grep -E "Transcript level|Intron chain|Matching intron" cmp.stats
```

### Per-bundle read-path diff

Reusable Python in memory file `project_rustle_parity.md` under "Read→transfrag
path instrumentation". Filters to a specific bundle's coord range, computes
shared / same-path / diff counts.

---

## 8. Memory file (~/.claude/.../memory/project_rustle_parity.md)

Full investigation history. Key sections:
- Session 1-3: 16 fixes establishing the baseline
- Session 4: graph/transfrag dump, root cause narrowed to sub-bundle decision
- Session 5: read→transfrag path instrumentation, characterized 4587 cases
- Session 6: dual-junction-set prototype, identified that "rescue everything"
  always regresses
- Session 7: architectural alignment with outer scope (this document)

**Always read the memory file before starting any new investigation** —
it documents what's been tried and why each approach failed/succeeded.

---

## 9. Default-on flags summary (post-session-7)

```
# Default-on (architectural):
#   - Outer junction stats in sub-bundles (initial seed + final filter)
#   - Skip sub apply_junction_filters_and_canonicalize
#   - Skip sub correct_junctions_with_map
#   - Skip sub snap_junctions_to_guides
#   - Skip sub good_junc_second_stage
#   - Skip sub apply_higherr_demotions / demote_runthrough / good_junc
#   - Skip sub repair_reads_after_junction_quality
#   - Inherit outer junction_redirect_map
#   - ST junction-strand color (RUSTLE_LEGACY_COLOR_BREAK to opt out)
#   - Exon-skip split (RUSTLE_EXON_SKIP_SPLIT_OFF to opt out)
#
# Recommended for highest Sn:
#   RUSTLE_ALT_SPLICE_RESCUE=1   (gives +0.4 Sn over default)
#
# Useful experimental (default-off, regress):
#   RUSTLE_GOOD_JUNC_STRICT, RUSTLE_LONGTRIM_STRICT_LEND,
#   RUSTLE_SHARED_SITE_WEAK_FILTER, RUSTLE_BADJUNC_TRUST_NODE_BOUNDARY,
#   RUSTLE_DUAL_JUNC_SET, RUSTLE_ALT_ACCEPTOR_RESCUE, RUSTLE_BUNDLE_GRAPH
```
