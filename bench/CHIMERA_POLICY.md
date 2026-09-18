# `--chimera-policy` — measured against `docs/PREREG_chimera_policy_2026-09-18.md` (md5 57f578ff0c4e2e4829d7c72bdfc5bb26)

**Verdict: the lever is well-defined and the chimera set is clean, but it changes no defensible number.
NOT ADOPTED as a reporting default.** Every apparent gain is disqualified by the pre-registered CP-5
denominator guard. Recorded so the question is closed with numbers instead of intuition.

Scripts `/mnt/linuxdisk/home/juanfraitu/chimera_policy/{arm_a,arm_b}.py`. `strict` on arm A reproduces
`soto_holdout.out` cell-for-cell — that is the correctness check on the re-implementation.

## The chimera set (frozen from the annotation before scoring)

| substrate | rule | count |
|---|---|---|
| human RefSeq | `description=` contains "readthrough" | 209 nodes; **1 of 27 NPIP members** (`PKD1P6-NPIPP1`); **0 of 19 TBC1D3 members** |
| CAT/GENCODE chr5/7/21 | name `A-B`, both halves are gene records | 137 genome-wide; **15 on chr5/7/21** |

## CP-1 — does chimera handling change the §6ks/§6kt conclusion? **ROBUST**

E1S beats E1 under all three policies, on both views:

| policy | view | bipartite F, E1 → E1S | pairwise prec, E1 → E1S | verdict |
|---|---|---|---|---|
| strict | universe | 0.831 → 0.881 | 0.815 → 1.000 | OK |
| strict | strict | 0.667 → 0.761 | 0.400 → 0.629 | OK |
| exclude | universe | 0.831 → 0.881 | 0.815 → 1.000 | OK |
| exclude | strict | 0.667 → 0.761 | 0.400 → 0.629 | OK |
| multi | universe | 0.831 → 0.881 | 0.815 → 1.000 | OK |
| multi | strict | 0.667 → 0.761 | 0.400 → 0.629 | OK |

The shipped `--min-shared-exon-frac 0.30` default is not conditional on chimera handling.

## CP-2 — arm A (held-out chr5/7/21): **NO-OP**, and structurally so

All twelve cells move by exactly 0.000. The reason is not that chimeras are harmless here — it is that
**none of them is in scope**: 0 of the 15 chr5/7/21 chimeras falls inside the 43-gene (E1) / 38-gene (E1S)
scored universe, and only 1 of the 15 (`C1QTNF3-AMACR`) is in Soto's 2,334-gene universe at all. The other
14 are ordinary protein-coding readthroughs (`PDCD6-AHRR`, `GIMAP1-GIMAP5`, `DUS4L-BCAP29`, …) that no
multi-copy family in the panel touches.

**Per the pre-registered CP-2 rule, the lever's scope is hereby restricted in writing to families that
actually contain a chimera.** Arm A cannot test it; it only shows the lever is inert where it does not apply.

## CP-3 / CP-5 — arm B (NPIP, TBC1D3): every gain is denominator shrinkage

| family | level | policy | genes | parts | R | P | **F** | pairwise sens | pairwise prec |
|---|---|---|---|---|---|---|---|---|---|
| NPIP | L1a/L1b | strict | 123 | 2 | 0.211 | 0.213 | **0.212** | 0.926 | 0.044 |
| NPIP | L1a/L1b | exclude | 111 | 2 | 0.225 | 0.227 | **0.226** | 0.923 | 0.050 |
| NPIP | L1a/L1b | multi | 123 | 2 | 0.217 | 0.219 | **0.218** | 0.932 | 0.047 |
| NPIP | L2 | strict | 91 | 2 | 0.286 | 0.289 | **0.287** | 0.926 | 0.081 |
| NPIP | L2 | exclude | 77 | 2 | 0.325 | 0.329 | **0.327** | 0.923 | 0.105 |
| NPIP | L2 | multi | 91 | 2 | 0.289 | 0.292 | **0.290** | 0.932 | 0.084 |
| NPIP | L3 | strict | 84 | 2 | 0.310 | 0.313 | **0.311** | 0.926 | 0.096 |
| NPIP | L3 | exclude | 74 | 2 | 0.338 | 0.342 | **0.340** | 0.923 | 0.114 |
| NPIP | L3 | multi | 84 | 2 | 0.303 | 0.307 | **0.305** | 0.929 | 0.092 |
| TBC1D3 | L1a/L1b | strict | 38 | 8 | 0.316 | 0.387 | **0.348** | 0.386 | 0.142 |
| TBC1D3 | L1a/L1b | exclude | 33 | 8 | 0.364 | 0.462 | **0.407** | 0.386 | 0.203 |
| TBC1D3 | L1a/L1b | multi | 38 | 8 | 0.333 | 0.406 | **0.366** | 0.414 | 0.159 |
| TBC1D3 | L2 | all three | 29 | 8 | 0.414 | 0.545 | **0.471** | 0.386 | 0.286 |
| TBC1D3 | L3 | all three | 19 | 10 | 0.526 | 1.000 | **0.690** | 0.263 | 1.000 |

**CP-3 (bar: |ΔF| ≥ 0.02 vs strict).**

| family | level | policy | ΔF | universe change | CP-3 | CP-5 guard |
|---|---|---|---|---|---|---|
| NPIP | L2 | exclude | **+0.040** | −15.4% | MOVED | ⚠ **> 10% ⇒ denominator shrinkage, not a gain** |
| NPIP | L3 | exclude | **+0.029** | −11.9% | MOVED | ⚠ **> 10% ⇒ denominator shrinkage, not a gain** |
| TBC1D3 | L1a/L1b | exclude | **+0.059** | −13.2% | MOVED | ⚠ **> 10% ⇒ denominator shrinkage, not a gain** |
| NPIP | L1a/L1b | exclude | +0.014 | −9.8% | no | (under the bar anyway) |
| NPIP | L1a/L1b/L2 | multi | +0.003…+0.006 | 0% | no | — |
| NPIP | L3 | multi | **−0.006** | 0% | no | `multi` can *hurt*: a second truth label the prediction cannot match |
| TBC1D3 | L1a/L1b | multi | +0.018 | 0% | no (0.018 < 0.020) | — |

**Every one of the four cells that cleared CP-3 was disqualified by CP-5.** `exclude` improves F only by
deleting 10–15% of the scored universe. `multi` never clears the bar, and at NPIP L3 it makes things
slightly worse — the honest cost of giving one record two truth labels.

## CP-4 — does `multi` open a DNA certificate? **NO — exactly as pre-declared**

Every (h_join, h_split] is empty at every level under every policy. h_split is −inf throughout: the truth
family is not even connected at the shipped cut, so no threshold can isolate it.

| family | level | h_join (all policies) | h_split | top boundary edges |
|---|---|---|---|---|
| NPIP | L1a | 1.0000 | −inf | CLN3@1.0000, EIF3CL@1.0000, PKD1P2@0.9994 |
| NPIP | L1b | 1.0000 | −inf | CLN3@1.0000, EIF3CL@1.0000, LOC100190986@1.0000 |
| NPIP | L2 | 1.0000 | −inf | CLN3@1.0000, LOC100190986@1.0000, LOC124907830@1.0000 |
| NPIP | L3 | 1.0000 | −inf | CLN3@1.0000, LOC100190986@1.0000, PKD1P5-LOC105376752@0.9994 |
| TBC1D3 | L1a | 1.0000 | −inf | NPEPPSP1@1.0000, TBC1D3P1-DHX40P1@0.9882, TBC1D29P@0.9614 |
| TBC1D3 | L2 | 1.0000 | −inf | USP6@1.0000, TBC1D29P@0.8119, LINC01743@0.2931 |
| TBC1D3 | L3 | 0.9614 | −inf | TBC1D29P@0.9614, USP6@0.9280 |

`exclude` changes exactly one visible boundary edge (NPIP L3: `PKD1P5-LOC105376752` → `LOC124907845`,
0.9994 → 0.9993) and leaves h_join at 1.0000. This is §6m1 reproduced through a different mechanism: the
binding constraint is **CLN3, EIF3CL and the LOC lncRNAs at identity 1.000**, none of which is a chimera.

## Reading

The lever does what it was asked to do — it makes "count the chimera" and "don't count it" directly
comparable, from a curated annotation-side label that never touches a prediction — and the answer it
returns is that chimera handling is **not** where the family definition is losing. The numbers that move
are the ones CP-5 was written to catch.

Kept as a **reporting option, off by default**. It is worth running whenever a claim is made about a family
that contains a chimera (NPIP does; TBC1D3 does not), precisely so that claim can be shown not to depend
on the policy.
