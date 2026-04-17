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
