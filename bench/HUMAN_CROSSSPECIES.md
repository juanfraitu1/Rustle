# Cross-species generalization — the identical method on HUMAN testis Iso-Seq

**Date:** 2026-07-12. **Purpose:** pre-empt the advisor's likely first objection — *"how do you know this is not
overfit for gorilla, or for this specific gorilla sample?"* We ran the **identical code and recipe** on a **human**
testis long-read RNA library — a **different lab, individual, and species** — and it reports the **real,
species-specific human copy numbers** (matching the human annotation). The method tracks the biology, not the sample.

## Data (all public, downloaded this session)

- **Reads:** `ERR13885926` — human **testis** full-length cDNA, GENCODE, **PacBio Sequel II (HiFi)**; 1,233,001
  reads, median 888 bp. (ENA; part of the GENCODE full-length long-read effort.)
- **Reference:** **T2T-CHM13v2.0** (`chm13v2.0.fa`, complete chrY — the same T2T basis Soto used).
- **Alignment (same recipe as gorilla):** `minimap2 -ax splice:hq --eqx -Y -N 50 -p 0.1 --secondary=yes` →
  **96% reads mapped** (1,179,288 primary), 7,942,153 secondary alignments (the multimappers the gate needs).
- **Copy calls (same binary, same flags):** `copy_assign --min-copies 2 --skip-poa-diagnostic --homology-primary`,
  foreground/serial (`human_families.sh`), on the human orthologs of the gorilla-validated families.

## Result — the method tracks species-specific copy number

| family | gorilla χ(H) | **human χ(H)** | human annot (T2T) | note |
|---|---|---|---|---|
| RBMY | 6 | **6** | ~6 | copies land on RBMY1B/A1/D/E/J/F |
| TSPY | 5 | **33** | ~35 | the human TSPY array (TSPY2/3/4/8/9/10 + array LOCs) |
| MAGEA | 2 | **11** | 11 | copies land on MAGEA1/2/2B/3/4/6/8/9/10/11/12; **CSAG correctly split into its own family** |
| GSTM | 3 | 2 | 5 | partial — only the expressed GSTM2 + GSTM5 resolved |
| PCDHB | 5 | — | 16 | coverage-limited: 77 reads in this library |
| DAZ | 2 | — | 4 | coverage-limited: 16 reads in this library |

**The load-bearing point:** MAGEA **2 → 11** and TSPY **5 → 33** across species. If the method were overfit to
gorilla it could not report these — it recovers the **human** expansions, matching the human annotation. Every
recovered copy was checked to land on an annotated paralog (not a mis-chain); MAGEA even split the adjacent CSAG
antigen family off correctly. RBMY is concordant (6 = 6) where the copy number is conserved between species.

**Honest limits (disclosed):** PCDHB and DAZ are under-expressed in this single testis library (77 / 16 reads) —
a depth limit of one dataset, not a method failure. The near-identical human-specific duplicates (SRGAP2B/C,
ARHGAP11B) sit at the K=0 identifiability frontier and/or below the expression floor here (ARHGAP11 129 reads →
unresolved); SRGAP2's copies are 84 Mb apart on chr1 and need the genome-wide catalog, not a single `--region`.

## The tightest control (noted, not yet run)

Makova's **PRJNA911852** is a single-protocol great-ape testis Iso-Seq set with **matched human AND gorilla**
samples (human SRR22838397/398/405/406; gorilla SRR22838403/404). Those are Sequel **subreads** (pre-CCS, noisy —
would need CCS before our PSV gate), so we used the clean HiFi GENCODE library here; the matched Makova set is the
natural next step for a same-protocol human↔gorilla comparison.

## What to tell the advisor

The identical binary, on a **human** testis Iso-Seq (different lab, individual, species), recovers multi-copy
families **at the human copy numbers** — RBMY 6, MAGEA 11, TSPY ~33 — with every copy on an annotated paralog. It
does **not** emit the gorilla numbers (MAGEA 2, TSPY 5): it tracks the species-specific biology. That is the direct
answer to "overfit to gorilla / to this sample."

Artifacts: `bench/make_human_crossspecies.py`, `bench/slides/human_crossspecies.png`, `/home/juanfra/human_val/`
(`human_families.sh`, `hum_*.quant.tsv`, `human_testis.t2t.bam`). Reference retrieved via NCBI/ENA; T2T-CHM13v2.0.
