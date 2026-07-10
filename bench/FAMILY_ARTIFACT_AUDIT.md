# Are we calling any more artifacts "families"? A catalog audit

**Date:** 2026-07-10. Current binary (`a190a77`). Ran `copy_assign` on **30 annotated multi-copy gene-family
loci** (2–12 paralogs, span < 900 kb, named families sampled from `GGO_genomic.gff`) plus the session panel,
foreground/serial, and classified every family that formed. The dense LOC/ZNF mega-clusters (>25 copies, >1.5 Mb)
were excluded — they time out and would over-form by repeat content regardless.

## What formed

Most named families are **silent in testis** (HOX, TAS2R, S100A, MMP, CLEC, SERPINA/B, KRT, WFDC, DEFB…): 0
reads, 0 families — correct. **16 families** formed across the loci that are expressed:

| class | n | families |
|---|---|---|
| **Clean** (one named gene family, copies span-disjoint) | 10 | GBP (4+2), PCDHB (5), APOBEC (2), APOL (3), TCEAL (3+2), GSTM (3), MAGEA (2), DAZ (2) |
| Real, LOC-named (annotation uses LOC symbols) | 2 | TSPY array (5), an r4 lncRNA pair |
| **Coverage-floor Containment** (known, flagged not pruned) | 2 | RFPL, r4 — `bench/CONTAINMENT_COVERAGE_FLOOR.md` |
| **NEW artifact: KRAB-ZNF domain bridge** | 1 | ZNF445 / ZKSCAN7 / ZNF197 |

## The one new artifact: a KRAB-ZNF domain bridge, and it is `--homology-primary`-only

`NC_086017.1:54339728-54640140` emitted a 3-copy family whose members are **three distinct annotated genes** —
ZNF445, ZKSCAN7, ZNF197. They are not a recent paralog family:

- Genomic spans align at **43–68% identity over ~1% coverage** — only the shared zinc-finger exon.
- Their **spliced exon-sum transcripts do not align at asm20 at all** (the 0.80 tier). The edge was formed by
  E_r's **sensitive tier (0.60 identity)**, which is meant to recover ancient paralogs but also bridges the
  zinc-finger domain shared across the whole KRAB-ZNF superfamily.

**Under the default E_c conflict oracle this cluster forms 0 families** (0 conflict edges — the three genes map
uniquely, so no read is ambiguous between them). The over-merge is therefore **specific to the opt-in
`--homology-primary` path**, and it comes from the E_r sensitive tier.

This is the exact domain-bridge trap the genome-wide family-definition work fought (the "238-blob" was KRAB-ZNF
domain bridges), and it was solved there by co-threading + community de-bridging on top of the raw E_r edges.
`copy_assign`'s inline `--homology-primary` branch (denovo_pipeline.rs:1256–1260) uses `homology_edges_all_reps`
followed by a **plain `gamma_quasi_clique_partition(0.20)`** — the raw E_r edges with **no de-bridging**. So it
recovers real ancient paralogs (GWFAM8 KRAB-ZNF, GSTM's globin problem) *and* the domain bridges, indiscriminately.

## Verdict

- **The default (E_c) catalog is clean of over-merges on this panel.** Every family it forms is a single named
  gene family with span-disjoint copies. The two flagged Containment cases (RFPL, r4) are the coverage floor from
  the previous turn, already reported-not-pruned.
- **`--homology-primary` has one artifact class: KRAB-ZNF (and by extension any large domain-sharing superfamily)
  domain bridges**, admitted by the E_r sensitive tier because a zinc-finger transcript is mostly zinc-finger
  domain. It is opt-in, so it does not affect the shipped default, but it is real whenever `--homology-primary`
  is used.

## Two minor notes (not artifacts, but worth recording)

- **Family fragmentation.** GBP (6 annotated) formed as 4 + 2 copies in two families; TCEAL (6) as 3 + 2. This is
  *under*-merging (a real family split), the opposite of the domain bridge — less harmful (no false copies) but
  it means a copy-number readout of "4" or "2" understates the true 6. Cause not investigated.
- The 30-locus sample is expression-limited: only ~6 of 30 named families are testis-expressed enough to form
  anything. A true genome-wide precision number needs the dense clusters (excluded here for cost) and is unmeasured.

## The fix, if wanted

The domain bridge is a **harmony gap** between `copy_assign --homology-primary` and the genome-wide E_r catalog:
the inline path skips the de-bridging the catalog path applies. The options, in order of scope:
1. Gate the E_r sensitive-tier edges in `copy_assign` on reciprocal **coverage** (the domain bridge is ~1%
   coverage; a real paralog pair is > 0.50) — cheap, threshold-in-the-existing-refine-gate, and it is exactly
   what `RefineParams` already carries but the sensitive tier apparently relaxes.
2. Port the co-threading / community de-bridging from the genome-wide catalog into the inline branch (DRY: call
   the same decomposition both places).

Option 1 is the minimal fix and the one to validate first: confirm the ZNF bridge is a low-coverage edge and that
raising the sensitive-tier coverage floor drops it without dropping GSTM / DAZ / the real KRAB-ZNF families.

Related: `bench/CONTAINMENT_COVERAGE_FLOOR.md`, `bench/YAG_VS_ISOCON.md`, `project_rna_family_homology_primary`,
`project_family_def_readconflict` (the domain-bridge history).
