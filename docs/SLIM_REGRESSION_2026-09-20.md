# Does the slimmed codebase still produce the recorded numbers?

Run 2026-09-20 at `3566682a`, after §6s0–§6s5 removed 3,884 lines of Rust and archived 251 scripts.
Every row was **re-executed today** against the current binaries; "expected" is the number already
written down in `REPRODUCE.md`, `bench/LAB_DATASET_BAKEOFF.md` or `docs/PIPELINE_STATE_TEST_2026-09-19.md`.

## Assembly — all four arms exact, and the head-to-head is unchanged

Scored with `gffcompare -r <ref>`; chains = matching intron chains.

| arm | mRNAs | chains | chain Sn/Pr | tx Sn/Pr | vs recorded |
|---|---|---|---|---|---|
| chr20 shallow (`human_testis`), polished | 658 | 336 | 7.8 / 51.6 | 7.4 / 51.2 | **exact** |
| A119b chr20, polished | 5,844 | 1,064 | 24.8 / 19.9 | 23.3 / 18.2 | **exact** |
| GGO `NC_073244.2`, polished | 4,063 | 1,575 | 28.4 / 38.9 | 26.6 / 38.8 | **exact** |
| GGO, k3 + junction-majority + polish | 5,298 | 1,688 | 30.4 / 31.9 | 28.5 / 31.9 | **exact** |

The competitors were re-scored from their stored GTFs against the same reference, so the comparison is
not carried over from a doc — it was recomputed:

| tool (GGO `NC_073244.2`) | mRNAs | chains | chain Sn/Pr | tx Sn/Pr | vs recorded |
|---|---|---|---|---|---|
| StringTie | 3,735 | 1,374 | 24.7 / 37.0 | 23.2 / 36.8 | **exact** |
| FLAIR | 6,169 | 1,393 | 25.1 / 23.6 | 23.5 / 22.6 | **exact** |
| isoseq collapse | 20,643 | 1,655 | 29.8 / 9.0 | 28.1 / 8.1 | **exact** |

⭐ The two standing claims hold on the slim tree: the shipped polish **beats StringTie and FLAIR on all
five metrics**, and k3+majority+polish **beats isoseq collapse outright** (1,688 chains at 31.9% precision
against isoseq's 1,655 at 9.0%, from a quarter of the transcripts).

## O1 guided — chr16, every number exact

RefSeq CHM13 chr16 gene+pseudogene bodies → `minimap2 -x asm20 -c --eqx -P` all-vs-all →
`mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`.

| | expected | got |
|---|---|---|
| graph | 380 nodes, 627 edges | **380 / 627** |
| pairs dropped, no exonic evidence | 10,637 | **10,637** |
| pairs dropped, shared-exon fraction < 0.60 | 965 | **965** |
| families (≥ 2), members, largest | 90, 271, 26 | **90, 271, 26** |
| size distribution | 56×2 19×3 6×4 3×5 2×6 2×7 1×11 1×26 | **identical** |

⭐ **NPIP: all 21 chr16 NPIP genes in ONE cluster (MCL0), none elsewhere** — NPIPA1/2/5/6/7/8/9,
NPIPB2–B15, NPIPB10P, NPIPB14P, with the same five unnamed LOC members (LOC100505915,
LOC124907807/808/834, LOC128966608). MCL1 is SMG1 + SMG1P1…P7; the PKD1P readthroughs
(PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P6-NPIPP1) cluster **separately** from NPIP, as recorded.

## O1 de novo — chr16, every number exact

A119b chr16 (1,787,427 records / 682,958 primary) → `copy_assign --assemble-only` + shipped polish →
de novo loci → all-vs-all → `mcl_families`, same flags. Annotation used **only to score**.

| | expected | got |
|---|---|---|
| transcripts | 9,629 | **9,629** |
| de novo loci | 2,550 | **2,550** |
| graph | 864 nodes / 1,836 edges | **864 / 1,836** |
| families (≥ 2) | 70 | **70** |
| NPIP genes covered | 21/21 | **21/21** |
| dominant cluster covers | 20/21 | **20/21** |
| NPIP-touching loci / clusters | 37 in 9 | **37 in 9** |
| median de novo loci per NPIP gene | 2 | **2** |
| the one gene outside the dominant cluster | NPIPB5 | **NPIPB5** |

## TBC1D3 — coherent, but no exact prior baseline for this run

chr17 guided, same recipe (430 nodes / 478 edges, 90 families, largest 12). **11 TBC1D3 copies in one
cluster** (TBC1D3, B, D, E, F, G, H, I, K, P1, P2) together with **USP6**, the progenitor TBC1D3 arose
from — the family is *unshattered*, which is what §6ks/§6kt claimed for `min_shared_exon_frac`. The
`TBC1D3P1-DHX40P1` readthrough clusters separately (with DHX40), the same boundary behaviour NPIP shows
with PKD1P; TBC1D3P3/P4 form their own pair. ⚠ `docs/PIPELINE_STATE_TEST_2026-09-19.md` tested chr16
only, so there is **no recorded chr17 figure to match against** — this row is a fresh measurement, not a
reproduction, and should not be quoted as "unchanged".

## One apparent regression that is not one

Re-running the **2026-09-07** human NPIP catalog (`soto_mcl/npip_hsa/cat.log`) gives 139 clusters / 424
members / largest 26, against that log's 173 / 505 / 27. The cause is in the log itself: the current
binary prints `min-shared-exon-frac=0.3: 761 pair(s) dropped`, a line the Sep 7 run does not have.
`min_shared_exon_frac` was promoted to the shipped default (0.30) in **`cdf1b4d0` (§6kt, a user
decision)** — *after* that log was written, and it is the change §6ks measured as bipartite F 0.83→0.88
with precision 0.82→1.00. The Sep 7 numbers are a superseded definition, not a target. No commit in
§6s0–§6s5 touches family-definition code.
