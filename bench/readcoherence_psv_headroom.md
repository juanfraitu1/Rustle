# Read-coherence × PSV — HEADROOM probe (go/no-go)

Cheap geometric+structural probe: do read-coherence's recall wins and PSV copy-resolution
act on the SAME loci (→ unified molecule-threaded graph worth building), or disjoint loci
(→ not worth the rebuild; keep them as two separate levers)?

## Inputs
- read-coherence recall set: transcripts tagged `source "read_coherence"` in rcg_*.gtf (25 chroms)
- multi-copy families (universe.tsv, n_copies>=2): **58** families, 556/614 member transcripts located in RefSeq
- family regimes: TANDEM_NEAR=27, DISPERSED=17, COLLAPSED=13, SINGLE_LOCUS=1

## Recall set
- read_coherence transcripts total: **23791**
- multi-exon (the recall lever): **23791**
- of those, FSM (intron chain == a RefSeq transcript = real recall): **1921** (8.1%)

## Tier 1 — are the recall wins even at multi-copy family loci?
- recall isoforms overlapping a multi-copy family locus: **88** / 23791 (**0.37%**)
- of those, FSM: **7** (of 1921 FSM = 0.36% of all read-coherence real recall)

## Tier 2 — of the family-locus hits, which regime (does PSV actually SPLIT?)

| regime | recall isoforms | FSM | PSV leverage |
|---|---|---|---|
| COLLAPSED | 35 | 4 | **YES — copies share one frame, PSVs split** |
| TANDEM_NEAR | 32 | 2 | no — distinct coords separate them |
| DISPERSED | 21 | 1 | no — distinct coords/contigs (RABL2 boundary) |

## Tier 3 — of the COLLAPSED hits, are the 'copies' REAL paralogs or domain-sharers?
(PSVs can only copy-resolve genuine paralogous copies. Two different genes that merely share a
protein domain/repeat and happen to overlap by annotation are NOT copies — PSVs are meaningless there.)

- COLLAPSED recall isoforms at CONFIRMED real paralog families (APOBEC3/RFPL/RABL2): **0**
- COLLAPSED recall isoforms at DOMAIN-SHARER false families (CDPF1/PPARA, CREB1/METTL21A, GCA/KCNH7, TTN/NCAPH2/MIEF1-LOC): **35**
- COLLAPSED recall isoforms at unclassified families: **0**

- read-coherence recall isoforms landing on ANY confirmed-real paralog family (any regime): **3**

## Verdict
Geometric PSV-resolvable headroom (recall isoforms at COLLAPSED loci) = **35** (0.147% of the recall set, 4 of them FSM).
TRUE headroom after removing domain-sharer false families = **0** (confirmed-real collapsed paralog recall isoforms = 0).

**GO/NO-GO: NO-GO** for a molecule-threaded graph justified by *combining* recall + copy-resolution.
Read-coherence recall (99.6% at single-copy loci) and PSV copy-resolution (needs collapsed real
paralogs) are disjoint in this data. The recall lever already ships additively (`--read-coherence`);
threading PSVs through the same graph adds copy-resolution to ~0 real recall isoforms.

Top families receiving recall isoforms:
  - LOC134756662 (COLLAPSED, 2 copies): 14 isoforms [DOMAIN-SHARER]
  - CREB1 (COLLAPSED, 2 copies): 7 isoforms [DOMAIN-SHARER]
  - LOC109024534 (TANDEM_NEAR, 2 copies): 6 isoforms
  - CDPF1 (COLLAPSED, 2 copies): 5 isoforms [DOMAIN-SHARER]
  - IGLL1 (DISPERSED, 2 copies): 5 isoforms
  - LOC101126655 (TANDEM_NEAR, 11 copies): 5 isoforms
  - GCA (COLLAPSED, 2 copies): 4 isoforms [DOMAIN-SHARER]
  - GGT1 (DISPERSED, 2 copies): 4 isoforms
  - LOC101136027 (TANDEM_NEAR, 2 copies): 4 isoforms
  - LOC129529456 (COLLAPSED, 2 copies): 4 isoforms [DOMAIN-SHARER]
  - CCDC188 (DISPERSED, 4 copies): 3 isoforms
  - LOC101129569 (TANDEM_NEAR, 10 copies): 3 isoforms
  - LOC101134642 (DISPERSED, 2 copies): 3 isoforms
  - LOC115931965 (TANDEM_NEAR, 2 copies): 3 isoforms
  - LOC129529434 (DISPERSED, 4 copies): 3 isoforms

![funnel](readcoherence_psv_headroom.png)

## What this does and does NOT settle
- **Settles:** the *PSV-unification* motivation. Read-coherence's recall win does not co-occur
  with PSV-resolvable copies, so building one graph to deliver BOTH on the same loci pays ~nothing.
- **Does NOT settle:** the molecule-threaded graph as a *pure recall* architecture (gate-not-kill,
  keep per-molecule evidence through paths). That stands or falls on recall alone — and the recall
  is *already* delivered additively by the shipped `--read-coherence` layer, so a rebuild would buy
  cleaner architecture / less noise, not new recall. It would NOT add copy-resolution.

## Honest caveats (which way they push)
1. **Geometric undercount (pushes headroom UP a little):** COLLAPSED is detected only when *annotated*
   copies overlap. Cross-mapping can pile a divergent/unannotated copy's reads onto a *single-copy*
   annotated locus, creating PSV signal the annotation hides. But prior copy_split real-data work
   (RABL2 / DAZ / co-located tests) showed that regime is rare and coverage-limited in GGO; even a
   3–5× undercount keeps the headroom sub-1% and the *real-paralog* count near zero.
2. **Small family universe (neutral):** 58 multi-copy families is the genome's reality, not a sampling
   gap. Multi-copy paralogs are a small slice of loci; read-coherence's recall is dominated by
   single-copy alt-splicing/novel-junction isoforms regardless of family-set size.
3. **Domain-sharer contamination (pushes headroom DOWN, decisively):** all 35 COLLAPSED hits land on
   families the Compara + contiguous-core validation already flagged as domain-sharing FALSE families
   (two different genes sharing a repeat, overlapping by annotation accident) — not paralog copies.
   PSVs cannot copy-resolve non-paralogs, so the TRUE headroom is 0, not 35.

## Reproduce
- `python3 bench/readcoherence_psv_headroom.py`  (reads /tmp/gw/rcg_*.gtf + ref_*.gff3 + universe.tsv)
- `python3 bench/readcoherence_psv_headroom_fig.py`  (renders the funnel figure)
- COLLAPSED-regime hit loci: `bench/headroom_loci.tsv`
