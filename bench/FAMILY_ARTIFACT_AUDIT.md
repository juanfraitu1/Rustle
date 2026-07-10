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

## The coverage fix was TESTED and it FAILS — reverted

The natural fix — **reciprocal coverage** (`min(qcov, tcov)` instead of coverage over the shorter transcript) —
was implemented and run on the panel. It correctly shrank the ZNF bridge (3 → 2 copies), but it **destroyed four
real families**: GSTM (3 → 0), MAGEA (2 → 0), PCDHB (5 → 0), TCEAL (5 → 0), and dropped GBP (6 → 4). Reverted.

The reason is fundamental. Reciprocal coverage penalises *transcript-length differences*, which are biologically
normal — isoforms, UTR extensions, and truncated de-novo models (DAZ2 is itself a 70% model). A real paralog
pair where one copy is longer than the other has the alignment cover the short copy fully but the long copy
partly, so `min(qcov, tcov)` falls below 0.50 and the edge is dropped. The coverage-over-shorter metric was
chosen deliberately to be robust to exactly this.

**The measured lesson: no pairwise coverage/identity threshold separates a domain bridge from a length-divergent
real paralog.** The ZNF bridge (shorter-coverage 0.55, identity 0.69) sits in the same pairwise cell as real
sensitive-tier paralogs. This is precisely why the genome-wide family-definition work solved domain bridges with
**graph structure** — co-threading + community de-bridging (weighted Louvain) — not a pairwise gate.

And that graph fix **cannot be ported to `copy_assign`'s per-region path**: the ZNF over-merge is a 3-node
complete triangle (ZNF445–ZKSCAN7–ZNF197, each pair domain-homologous), which any γ-quasi-clique accepts as a
family. De-bridging works only with the *genome-wide* graph, where the whole KRAB-ZNF superfamily blob is visible
and Louvain can cut it into sub-communities. A single 300 kb window with three ZNF genes has no larger structure
to exploit.

## Disposition

- **Do not ship a coverage/identity threshold.** It breaks real length-divergent families (measured) and is the
  arbitrary threshold the advisor rejects.
- The KRAB-ZNF domain bridge is an **inherent limitation of the per-region `--homology-primary` mode**: it
  inherits raw E_r edges with no genome-wide context to de-bridge. It is opt-in and does not affect the default
  E_c catalog (which forms 0 families there).
- For a clean real-family catalog, the **genome-wide de-bridged path** (`detect_homology_catalog_genome_wide` +
  co-threading community detection) is the correct tool — `copy_assign --homology-primary` is a per-region
  assignment mode, not a family-definition catalog.

The one honest cheap improvement already in place is that the default oracle is clean; the domain bridge only
surfaces when a user opts into per-region E_r membership on a domain-sharing superfamily.

Related: `bench/CONTAINMENT_COVERAGE_FLOOR.md`, `bench/YAG_VS_ISOCON.md`, `project_rna_family_homology_primary`,
`project_family_def_readconflict` (the domain-bridge history).
