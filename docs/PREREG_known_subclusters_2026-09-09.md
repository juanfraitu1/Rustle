# PREREG — do our families recover the KNOWN subcluster structure of NPIP and TBC1D3? (advisor, 2026-09-09)

**Written after reading both papers and BEFORE looking at any of our clustering output for this question.**

## The external truth, quoted

### NPIP — Dishuck et al. 2025 (`npip_dishuck_nihpp-2025.02.04.636496v1.md`)
- Two major subfamilies **NPIPA** and **NPIPB**. ⭐ **"IGC occurs within but not between the two major
  subfamilies (NPIPB copies undergo IGC only with NPIPB but not NPIPA loci)"** (l. 239–240) ⟹ the A/B split is
  a REAL sequence divergence that gene conversion does not erase.
- Nested groups inside B: **B3–B5** and **B6–B9** ("at least one paralog from each of these larger
  subfamilies is always present", l. 226–227); **B12/B13** are one clade because there is insufficient genetic
  distance to separate them (l. 188).
- A human-specific clade **B3, B4, B5, B11, B12, B13** shares a final-coding-exon change (l. 416).
- **NPIPA1** is the ancestral paralog and acts only as an IGC donor, never acceptor (l. 241, 496).

### TBC1D3 — Guitart et al. (`tbc1d3_guitart_subclusters.md`)
- **cluster 1** and **cluster 2**, holding the majority of paralogs, **1.35 Mbp apart** in apes (l. 88, 135),
  plus older **orphan** copies distributed along chr17 (l. 234, 241).
- ⭐⭐ **"complete lineage-specific stratification … into distinct clades for human, Pan, gorilla …"**, taken as
  evidence of **recurrent duplication or gene conversion of all gene family copies in each lineage**
  (l. 148–151). ⟹ within human, the copies are ONE clade; **clusters 1 and 2 are POSITIONAL, not phylogenetic**.
- The 43 bp ORF deletion is shared by **all** cluster-1 and cluster-2 copies and absent from orphans (l. 241).

## ⭐ The two families therefore predict OPPOSITE outcomes — that is what makes this a test
| | known structure | is it a SEQUENCE split? | what our method SHOULD do |
|---|---|---|---|
| **NPIP A vs B** | two subfamilies, no IGC between them | **yes** | **separate them** — a homology method that cannot is missing real divergence |
| **TBC1D3 cluster 1 vs 2** | two positional clusters, IGC homogenises within lineage | **no** | **NOT separate them** — recovering a positional split from sequence alone would be suspicious |

## Predictions
| # | prediction | refuted by |
|---|---|---|
| **P1** | Somewhere on the inflation ladder, NPIP splits into groups that align with **A vs B** | no inflation separates A from B at better than chance |
| **P2** | TBC1D3 clusters 1 and 2 are **NOT** recovered as the first cut; the ladder's cuts cross the positional boundary | a cut cleanly reproduces cluster 1 vs cluster 2 |
| **P3** | TBC1D3's observed I = 5.0 cut (**{D,K} vs the other seven**) is a **sister-pair**, not a cluster boundary — D and K sit 11.6 kb apart and are both 12,610 bp | D and K fall in different positional clusters |
| **P4** | Pairwise identity among the 9 human TBC1D3 paralogues is **uniformly high** with no bimodality at the cluster boundary, consistent with IGC homogenisation | identity is clearly bimodal by cluster |
| **P5** ⚠ | For NPIP, the finer groups (**B3–B5**, **B6–B9**) require *more* inflation than the A/B split, if they appear at all | a finer group separates before A/B |

## ⛔ Prior result that must be reported either way
§6ew (2026-09-05, rows 693–695) already found **MCL keeps all 20 NPIP loci together to I = 4.0** and that
"identity/coverage do NOT separate" A from B. **If P1 fails, it is a REPLICATION of that negative, not a new
finding**, and the honest statement is that our method does not resolve the A/B subfamily split.

## Rules held
⚠ Gene symbols are used **only to read the result out**, never to build it (`project_gene_naming_trap`).
⚠ Human and gorilla are never pooled — NPIP is scored on human CHM13 chr16 and TBC1D3 on human CHM13 chr17.
⚠ Positional cluster membership for TBC1D3 is assigned from coordinates BEFORE seeing any clustering output.
