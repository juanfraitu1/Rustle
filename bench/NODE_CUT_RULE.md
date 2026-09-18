# The NODE CUT rule — measured against `docs/PREREG_node_cut_2026-09-18.md` (md5 2af3393070d6c2ded8db3cd89c0a6dc3)

**Verdict: SAFE BUT HARMFUL. NOT ADOPTED.** Cutting a chimeric record at its parent boundary loses no
family member (NC-1 passed everywhere), but it makes every other number worse, and it fails the deciding
NC-2 test at every level. The mechanism is identified and is a general constraint on node rules, not a
detail of this one.

Scripts `/mnt/linuxdisk/home/juanfraitu/node_cut/{cut,score,fresh}.py`; data `cuts.tsv`, `records_cut.tsv`.

## Correctness check

The cut re-emits records from the PAFs `family_cert` already stored — **no new alignment**. Of the 8,669
baseline records, the **7,444 that touch no cut node are byte-identical** after re-emission (0 differing,
0 added, 0 lost). Only the 1,225 → 1,191 rows involving a chimera changed.

## The cut itself (NC-5)

| | |
|---|---|
| readthrough primary nodes | 209 |
| **cut** (labels form exactly two contiguous parent runs) | **195** |
| `uncut`, left whole and reported | **14** |
| node set | 58,563 → 58,758 (+390 pieces, −195 originals) |

Higher than the 173/181 known before the run because the pre-registered complement rule (§1.2) handles the
records where only one half is a gene record — which is what makes `PKD1P6-NPIPP1`, `PKD1P4-NPIPA8` and
`PKD1P5-LOC105376752` cuttable at all.

## NC-1 — safety: **PASSED**

Pairwise sensitivity is unchanged at **every** level, for both families, under both truth mappings:
NPIP 0.926 → 0.926; TBC1D3 0.386 → 0.386 (L1a/L1b/L2), 0.263 → 0.263 (L3). No member is lost.
The pre-registered majority-exon-bp mapping and the name-matched alternative picked the **same piece in
every case**, so that ambiguity never bound.

## DEV arm — every other number moves the wrong way

| family | level | arm | universe | R | P | **F** | sens | prec |
|---|---|---|---|---|---|---|---|---|
| NPIP | L1a/L1b | baseline | 123 | 0.211 | 0.213 | **0.212** | 0.926 | 0.044 |
| NPIP | L1a/L1b | cut | 141 | 0.184 | 0.186 | **0.185** | 0.926 | 0.033 |
| NPIP | L2 | baseline | 91 | 0.286 | 0.289 | **0.287** | 0.926 | 0.081 |
| NPIP | L2 | cut | 109 | 0.239 | 0.241 | **0.240** | 0.926 | 0.056 |
| NPIP | L3 | baseline | 84 | 0.310 | 0.313 | **0.311** | 0.926 | 0.096 |
| NPIP | L3 | cut | 98 | 0.265 | 0.268 | **0.267** | 0.926 | 0.070 |
| TBC1D3 | L1a/L1b | baseline | 38 | 0.316 | 0.387 | **0.348** | 0.386 | 0.142 |
| TBC1D3 | L1a/L1b | cut | 34 | 0.353 | 0.444 | **0.393** | 0.386 | 0.188 |
| TBC1D3 | L2 | baseline | 29 | 0.414 | 0.545 | **0.471** | 0.386 | 0.286 |
| TBC1D3 | L2 | cut | 30 | 0.400 | 0.522 | **0.453** | 0.386 | 0.261 |
| TBC1D3 | L3 | baseline | 19 | 0.526 | 1.000 | **0.690** | 0.263 | 1.000 |
| TBC1D3 | L3 | cut | 20 | 0.500 | 0.909 | **0.645** | 0.263 | 0.818 |

NPIP: universe **+14.6% to +19.8%**, F **−0.027 to −0.048**, precision down at every level. The single
improving cell (TBC1D3 L1a/L1b, ΔF +0.046) comes with a **−10.5% universe** and is reported as such under
the NC-4 guard, not as a gain.

## NC-2 — the deciding test on the FRESH arm: **FAILED at every level**

Families named in the prereg before any score. Five were scorable in this graph (ID_305, ID_380, ID_453,
ID_454, ID_480 — ID_121 has only one node here).

| level | arm | universe | R | P | F | sens | prec |
|---|---|---|---|---|---|---|---|
| L1a/L1b | baseline | 28 | 0.250 | 0.304 | 0.275 | 0.667 | 0.074 |
| L1a/L1b | **cut** | **288** | 0.024 | 0.025 | 0.025 | 0.954 | **0.004** |
| L2 | baseline | 28 | 0.250 | 0.304 | 0.275 | 0.667 | 0.074 |
| L2 | **cut** | **224** | 0.031 | 0.032 | 0.032 | 0.941 | **0.005** |
| L3 | baseline | 24 | 0.292 | 0.368 | 0.326 | 0.632 | 0.100 |
| L3 | **cut** | **202** | 0.035 | 0.036 | 0.035 | 0.935 | **0.005** |

Bar was precision **+0.05 with no recall loss**. Measured: precision **−0.070 to −0.095**, universe
**+700% to +929%**. The apparent sensitivity rise (+0.275 to +0.304) is the artefact of a component that
swallowed everything — exactly the §6m0 trap about reading a merge as recall.

## NC-3 — certificates: **NO, as pre-declared**

Every (h_join, h_split] stays empty at every level, h_split = −inf throughout. NPIP h_join stays 1.0000.
**TBC1D3 L3 gets worse**: h_join 0.9614 → 0.9916.

## Why it fails — two mechanisms, both general

**1. Cutting does not separate; it doubles.** After the cut, **7 of 8** chimeras have **both** halves inside
NPIP's L2 component (BOLA2-SMG1P6, PDXDC2P-NPIPB14P, PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P6-NPIPP1,
SLX1A-SULT1A3, SLX1B-SULT1A4); at L3, 6 of 6. The "PKD1 half" is itself homologous to the NPIP-family
segmental duplication, so it does not leave. One outsider becomes two.

**2. A shorter node is an easier node — pieces become hubs.** Recomputing coverage against the piece's own
length shrinks the denominator, so a piece clears the coverage floor against more targets than the whole
record did:

| piece | exon bp | degree at L2 | whole record's degree |
|---|---|---|---|
| `PDXDC2P-NPIPB14P#1` | 2,626 | **50** | 49 |
| `PDXDC2P-NPIPB14P#2` | 1,568 | 33 | 49 |
| `PKD1P4-NPIPA8#1` | 1,223 | 36 | 39 |
| `PKD1P3-NPIPA1#2` | 2,072 | 36 | 38 |

19 piece nodes carry 282 degree between them; L2 edges rise 1,337 → 1,417. **The cut adds edges.**

This is the `min_cov_longer` invariant (§6cr) reappearing: a coverage floor must be denominated on the
LONGER sequence, because a short fragment covers most of itself and almost none of a real gene. Any node
rule that SHRINKS a node inherits this problem.

## What survives

- The **theory answer stands**: a record is in two families because its pieces are, and cutting keeps
  levels as node partitions, so T1, T1′ and T2 are untouched. The rule is well-defined and computable
  with no new alignment.
- The **applicability finding stands**: 195 of 209 readthrough records have a single clean cut point.
- What does not stand is the expectation that this improves anything. It does not, on either arm.

**Constraint for the next node rule:** it must not shrink a node without re-denominating coverage on the
original length, and it must be checked for whether both pieces remain inside the family's component —
that check alone would have predicted this result before any scoring.
