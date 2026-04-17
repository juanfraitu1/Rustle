# Remaining Rustle vs StringTie divergence on GGO_19

## Setup

- BAM: `GGO_19.bam` (chr19 PacBio IsoSeq, gorilla)
- StringTie: `stringtie -L`
- Rustle (baseline, no `--vg`): `rustle -L`
- Compare with `gffcompare -r stringtie.gtf rustle.gtf`

## Headline

- **1472 of 1839** StringTie transcripts recovered as exact `=` match (80.0%)
- **295 StringTie transcripts missed** (16.0%)
- Locus-level: **556/583** StringTie loci hit (95.4%); 11 fully missed loci

## Breakdown of the 295 missed transcripts

| Category | Count | Description |
|----------|-------|-------------|
| **Alt-splice in middle** | **118** | differs from a matched sibling in interior introns (exon skip, alt donor/acceptor) |
| Not overlapping any matched sibling | 66 | StringTie isoforms with no matched sibling at the same gene — includes the 11 fully-missed loci |
| Alternative TTS | 52 | same entire 5' intron chain as a matched sibling, different last intron (alt-polyA) |
| Alternative TSS | 36 | same 3' chain, different first intron (alt-promoter) |
| Single-exon | 23 | single-exon minor isoforms |

## Concrete examples

### STRG.1 (chr19:15395227-15569117, − strand)

StringTie emits 5 isoforms; Rustle emits 3. Miss breakdown:

| Ref tx | Exons | cov | Last exon | Match |
|--------|-------|-----|-----------|-------|
| STRG.1.1 | 9 | 17 | 15522277–15522602 | **= RSTL.285.1** |
| STRG.1.2 | 9 | 15 | 15568933–15569117 | **= RSTL.285.2** |
| STRG.1.3 | 9 |  3 | 15517453–15517765 | **missed** |
| STRG.1.4 | 8 |  4 | 15522277–15522602 | **= RSTL.285.3** |
| STRG.1.5 | 9 |  2 | 15536650–15536882 | **missed** |

Pattern: 4 alternative TTS variants (STRG.1.1, .2, .3, .5 — same 8-intron prefix,
different last exon); Rustle emits 2 of 4.

### STRG.25 (chr19:16910423-16959006, − strand)

StringTie emits 8 isoforms; Rustle emits 6. Miss: STRG.25.6 (cov=1) and STRG.25.8
(cov=2). Both are specific TSS + TTS combinations absent from Rustle's output even
though the constituent splice choices appear separately in other matched isoforms.

## Interpretation

Rustle's gap is **not at the locus-discovery level** (locus sensitivity is 95.4%).
It's at the **minor-isoform level within multi-isoform genes**:

1. **Max-flow decomposition truncates at the tail**. The heaviest few paths are
   extracted; low-abundance (cov=1–3) paths that share most edges with the major
   paths get their flow absorbed into the majors and never emit.

2. **Alternative last-exon combinatorics are not fully enumerated**. If a gene has
   4 distinct last-exon choices observed in reads, Rustle often emits 2 of them
   (the heaviest) and collapses the rest into the dominant TTS variant.

3. **Middle alt-splice events drop out when the minor path's flow is < 1**.
   StringTie has path-tracking heuristics (`keeptrf`, `checktrf`) that preserve
   seed candidates even after major extraction; Rustle ports those mechanisms
   but still loses 118 middle-splice minor isoforms.

## What the 11 fully-missed loci look like

Sample (see `/tmp/missed_by_rustle.txt`):

- Many are 1-2 exon minor genes with cov < 3 — below Rustle's default keeptrf/seed
  thresholds but above StringTie's.

## Mitigation options

| Approach | Estimated gain | Risk |
|----------|---------------|------|
| Lower cov-based single-exon threshold | ~10–20 tx | precision hit |
| Preserve seed abundance during flow decomposition (already applied for `trflong_seed`) | already in | — |
| Second-pass extraction on residual flow after major extraction | ~50–80 tx | likely novel loci cost |
| Alt-TSS/TTS enumeration at terminal splice sites | ~88 tx (52+36) | duplicate false positives possible |

Current baseline (1472 matches, 80.0% sensitivity) is already close to the
architectural ceiling imposed by the flow-decomposition paradigm shared with
StringTie. The remaining 295 is a long-tail problem — most individual paths
have cov=1–3 and require per-locus targeted fixes rather than a single lever.

---

## Updated divergence after the new-rep weak fix

After `a72a352` (preserve trflong_seed eligibility when tf replaces a keeptrf rep):

| Metric | Baseline | After fix |
|--------|----------|-----------|
| Matching transcripts | 1472 | **1556** (+84) |
| Sensitivity | 80.0% | 84.6% |
| Missed StringTie tx | 295 | 231 |
| Missed loci (full) | 11 | 11 |

### Re-categorization of the remaining 231 missed transcripts

| Category | Before fix | After fix | Δ |
|----------|-----------|-----------|---|
| Alt-splice middle | 118 | 113 | −5 |
| Alt-TTS | 52 | 42 | −10 |
| Alt-TSS | 36 | 25 | −11 |
| Not overlapping any matched sibling | 66 | 28 | −38 |
| Single-exon | 23 | 23 | 0 |

The biggest absolute drop is in "not overlapping" (66→28, −38): these are loci
where our earlier assembly produced nothing at all, and the new-rep fix
preserved representatives that got the locus at least partially assembled.

## Second-pass analysis on remaining misses — different patterns

Continuing the same diagnostic (pick miss → mini-BAM → trace → find
divergence point) on several remaining misses:

### STRG.271.4 (altMid, cov=1, 9 exons)

Rustle recovers 4 of 5 STRG.271 isoforms as `=` match (was 3 of 5
pre-fix — the fix helped here). The remaining miss (STRG.271.4) is an
exon-skipping isoform; its specific exon-skip pattern is carried by only a
handful of reads and its internal junctions don't form a distinct
keeptrf seed — they get absorbed by the parent isoform's chain.

This is a **genuine long-tail case**: the flow decomposition can't
enumerate every exon-combination for cov=1 seeds without producing many
false positives.

### STRG.15.1 vs STRG.15.2 (altTSS, cov=2 each, differ only in first exon end)

- STRG.15.1: first exon 16598646–16599149 (length 503)
- STRG.15.2: first exon 16598646–16599242 (length 596)

Same TSS, different first-intron DONOR (96 bp apart). Rustle collapses
these alternative donor sites in its graph — node segmentation treats
16599149 and 16599242 as one splice donor, so only one isoform emerges.

This is a **graph-structure (node granularity)** issue, not a flow issue.
Fixable by enabling finer splice-site discrimination at graph construction
time; orthogonal to the seed-eligibility bug just fixed.

### STRG.110.1 (cov=28, 6 exons, fully missed — no sibling matched either)

Reads spanning both the upstream gene (22150424–22152352) AND STRG.110
(22152353–22155355) cause Rustle's bundle builder to merge them into a
single 22149137–22155355 bundle. Pre-filter extraction produces 7
candidate transcripts spanning the full region; the **read-chain witness
filter** then kills 2 whose consecutive intron pairs aren't supported
by any single read, and pairwise/junction filters kill the rest. Final
output retains only transcripts at the upstream gene (22150424–22152352).

This is a **bundle-splitting** issue: when two adjacent genes have a
~1 bp gap between them and read-through transcription produces
bridging reads, Rustle merges them. StringTie's bundle-splitting
heuristics are more aggressive here.

## Summary of remaining divergence classes

| Pattern | Count (est.) | Fix type |
|---------|--------------|----------|
| Systemic new-rep weak bug | fixed: −84 | done |
| Graph-node segmentation (alt donor/acceptor <100 bp apart) | ~30–50 | graph construction tuning |
| Bundle merging across adjacent genes | ~10–20 | bundle splitter refinement |
| Long-tail minor isoforms (cov=1-2, exon-skip combinations) | ~130–180 | per-seed cost models, not a single lever |

The one systemic bug in the remaining list is the graph-node segmentation
for close alternative splice sites. The bundle-splitting one is
locus-specific. The long tail is inherent to flow decomposition.
