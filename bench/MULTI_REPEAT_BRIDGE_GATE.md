# Multi-repeat-bridge conjunction gate — a 3rd VG-native family-refinement gate

**Question.** Does the CONJUNCTION gate *(repeat-bridged **AND** no full shared exon)* remove the
DISCONNECTED repeat-bridge RNA-level false-positive class — the dominant residual over-merge
(34 / 83 E_p-impure blocks, 41 %) where the pairwise E_r oracle merged loci sharing **no full exon**
at id ≥ 0.70 (cross-component best-exon-id median 0.586 = a sub-exonic Alu/repeat bridge, not real
transcript homology) — while **sparing GSTM2 / MAGE** (they share a full exon) and real paralogs, and
**holding the named diploid-DNA oracle recall R_oracle = 50/57**? Or does it bleed into real families
the way plain exon-containment did (`EXON_CONTAINMENT_CRITERION.md` FALSIFIED cutting on exon-coverage:
it lost real single-gene paralogs 1:1 on fragmentation / UTR noise)?

**Answer: YES — it removes the repeat-bridged half of the class at zero recall cost on the named gold,
spares GSTM2/MAGE structurally, and does NOT bleed into real single-gene paralogs — but it is NARROW
and NOT strictly zero-collateral catalog-wide.** Verdict + recommendation below.

Artifacts (deterministic, `PYTHONHASHSEED=0`):
- `bench/multi_repeat_bridge_gate.py` — characterize (`char`) + gate (`gate`) stages.
- `bench/multi_repeat_bridge_gate.tsv` — 63-family characterize table (md5 `5e40411a`, byte-identical across runs).
- `bench/multi_repeat_bridge_gate.json` — gate sweep, two scopes, provenance, verdict (md5 `043b22d7`, byte-identical across runs).

Reused (not re-derived): `recombination_bridge_detector.py` (full-exon locus graph + exon-match tensor +
DISCONNECTED test + `cross_component_maxid`), `vg_repeat_catalog.py`/`.tsv` (canonical-minimizer node sets +
global multiplicity + cyclic flag), `family_er_pr.load_meta`, `family_level_pr_current.py` (eval ctx +
diploid-DNA oracle + relabel), `graph_def_refine_sweep.py` (block eval + multi-copy filter).

---

## The gate

```
CUT family  iff
   (NO full shared exon : ≥2 full-exon components  AND  cross-component best-exon-id < 0.70)
   AND
   (REPEAT-bridged      : ≥ C distinct cross-component shared VG minimizer nodes, each with
                          global multiplicity ≥ T)
   AND
   (same-gene guard     : never separate loci of the same annotated gene)
On CUT: replace the family by its full-exon components; <2-locus components drop via the multi-copy filter.
```

The full-exon components are exactly what remains after cutting the repeat-bridged cross-component edges.
`T` (multiplicity threshold) and `C` ("several repeats" count) are swept.

### Why the conjunction, not either half alone
- **Multiplicity alone** (the shipped `min_shared_mult ≥ 20` gate) can't be lowered without hitting the
  controls: GSTM2's internal Alu node has multiplicity **426** and MAGE 8–10, so a mult-only cut at any
  low `T` would shred GSTM2. The **no-full-shared-exon conjunct** blocks this — GSTM2/MAGE form a **single**
  full-exon component, so there is **no cross-component pair** and the gate structurally cannot fire.
- **Exon-containment alone** was FALSIFIED (`EXON_CONTAINMENT_CRITERION.md`) — it cut real single-gene
  paralogs that fragment (UTR noise, exon-skip). The **repeat conjunct + same-gene guard** spare those.

---

## Class signature (characterize stage, 63 families)

| class | n | shares full exon | med x-comp exon-id | med rep-bridges m5 / m8 | med max-cross-mult | med max-internal-mult | cyclic bridge |
|---|---|---|---|---|---|---|---|
| **DISCONNECTED_FP** | 35 | **0 %** | **0.586** (== the MD figure) | **3 / 1** | **13 (sub-20)** | 17 | 60 % |
| **GSTM2** (fid 9,13) | 2 | **100 %** | — (single comp) | 0 / 0 | 0 | **426** | 0 % |
| **MAGE** (fid 546,548,553) | 3 | **100 %** | — | 0 / 0 | 0 | 10 | 0 % |
| **REAL_PARALOG** (controls) | 23 | **82.6 %** | 0.547 | 0 / 0 | 0 | 31 | 8.7 % |

The signature separates cleanly: DISCONNECTED_FP has `frac_no_full_shared_exon = 1.0` with sub-20
cross-component repeat nodes; GSTM2/MAGE are connected single components (gate can't fire) but GSTM2's
internal mult-426 Alu proves a mult-alone gate would over-cut it.

**Provenance fix (had to reconcile, not just reuse).** The shipped `recombination_bridge_detector.tsv`
(pre-split 606-fam) family_ids no longer index the current 605-fam post-split `family_rna_refine.tsv`
(only 13/34 fids overlap). DISCONNECTED was **re-derived on the current catalog** by reusing
`R.analyze_family`, reproducing the 34-block count. Roster = 35 = 34 re-derived DISCONNECTED
(33 genuine multi-gene over-merges + fid 554 = same-gene `LOC101152062×3` fragmented real, spared by the
guard) + the FOXO1 hub fid 331 forced in. GSTM2 = fid 9,13; MAGE = 546/548/553 (shifted −1 from stale
547/549/554).

---

## Gate result vs the SHIPPED eval machinery

Baseline (605 multi-copy families): impure **79**, P_Ep **0.8694**, **R_oracle 50/57 = 0.8772**,
kept_REAL 423/460, kept_GENUINE 1224, fp_multifam 4.

Two scopes: **ROSTER** (cut only the pre-identified DISCONNECTED-FP roster = the task's literal target)
and **CATALOG** (the same predicate deployed blindly on all 605 families = stress test).

### Sweep (T × C), ROSTER scope

| T , C | families cut | **DISCONNECTED removed** | impure (Δ) | **P_Ep (Δ)** | **R_oracle** | recov_mc | kept_REAL Δ | GSTM2 / MAGE / real-ctrl cut |
|---|---|---|---|---|---|---|---|---|
| **5 , 1** | 24 | **24 / 35** | 57 (−22) | **0.9034 (+0.034)** | **50/57 held** | 50 | −3 | 0 / 0 / 0 |
| 5 , 2 | 20 | 20 / 35 | 60 (−19) | 0.8985 (+0.029) | 50/57 held | 50 | −3 | 0 / 0 / 0 |
| 8 , 1 | 17 | 17 / 35 | 63 (−16) | 0.8939 (+0.025) | 50/57 held | 50 | −2 | 0 / 0 / 0 |
| **8 , 2** | 15 | 15 / 35 | 65 (−14) | 0.8909 (+0.022) | 50/57 held | 50 | −1 | 0 / 0 / 0 |
| 12 , 1 | 16 | 16 / 35 | 64 (−15) | 0.8924 (+0.023) | 50/57 held | 50 | −2 | 0 / 0 / 0 |
| 12 , 2 | 13 | 13 / 35 | 67 (−12) | 0.8880 (+0.019) | 50/57 held | 50 | −1 | 0 / 0 / 0 |

**Invariants across the ENTIRE grid (6 points × both scopes):** R_oracle held (recov_mc never < 50);
GSTM2/MAGE cut = 0; 23 real-paralog controls cut = 0; P_Ep improves. Confirmed deterministic.

### Best (T,C) and the E_p delta
- **Max removal = (T=5, C=1): 24/35 DISCONNECTED removed, P_Ep +0.034 (0.8694 → 0.9034), R_oracle 50/57 held,
  GSTM2/MAGE/real-ctrl 0, kept_REAL −3.** This is C=1 (a *single* mult≥5 cross-component node).
- **Multi-repeat operating point = (T=8, C=2): 15/35, P_Ep +0.022, kept_REAL −1** — quote this when
  invoking the "several repeats" name (C≥2). The 24/35 headline is the C=1 upper bound.

The 24 removed at (5,1), **named** — all functionally-unrelated gene pairs (genuine over-merges):
PDIA3+PRKAB2 (7), the 20-gene Alu hub **fam17** (17), FGD3+HIVEP1 (24), KYAT3+RBMX (30), DCP1B-hub (57),
DLD+EPPIN-hub (64), AKR1C3+IL20RB (165), B3GAT2+SMAP1 (182), CEP57+MTMR2 (258), ACAD10+SBK3 (278),
SPSB2+TPI1 (284), RHBDF1-hub (310), ATP6V1C2+PDIA6 (322), ACP1+ALKAL2 (324), the ANKRD18/FOXO1 hub (331),
ASB6+NTMT1 (350), EIF1AX+LOC (372), HERC2-hub (403), APOL-hub (419), HDGFL3+TM6SF1 (424), GSPT2+SNX29 (451),
MTREX+PLPP1 (485), ECSIT+EXOG+RPL18A (492), CAND2+RPL32 (586).

The 11 abstained at (5,1): the **17 low-mult-glue** cases (max-cross-mult ≤ 8: e.g. MOV10+RHOC (20),
ESRRA+TPTE2 (268), JMJD8+STUB1 (454), NKIRAS1+UBE2E1 (587)) — correctly ABSTAINED, they are not high-copy
repeats (candidate readthrough / micro-homology / low-mult paralog); plus the **3 same-gene-guard-blocked**
(e.g. SERF2+SERINC4+TTPAL (128) rep_m8=69 but one gene spans components → recall-safe spare) and fid 554
(same-gene `LOC101152062` fragmented real).

---

## GSTM2 / MAGE spared — and by the correct mechanism

Cut **0** at every (T,C) in **both** scopes. Mechanism: they share a full exon → single full-exon
component → **no cross-component pair exists** → the conjunction structurally cannot fire. GSTM2's
internal Alu node has multiplicity 426; a multiplicity-alone gate lowered to catch the sub-20 DISCONNECTED
bridges *would* over-cut it. The no-full-shared-exon conjunct is exactly what saves it.

## R_oracle 50/57 held — but it is insensitive by construction (honest caveat)

R_oracle (the named diploid-DNA gold) held at every grid point in both scopes: recov_mc = 50, and the
recovered-gene **set** is identical at every cut point (not a masked relabel swap). But the gene-level DNA
oracle measures per-**gene** recovery, and the **same-gene guard structurally preserves per-gene recovery**
— so R_oracle *cannot* detect this gate's only failure mode (cross-gene paralog splitting). The guard
protects exactly what the oracle measures. **The operative recall probe of the real failure mode is
therefore the conservative cDNA-pair truth `kept_REAL`, which does move — not R_oracle.** We report both;
we do not lead with "named gold held" as the sole safety proof.

## kept_REAL cost + provenance (the honest bleed)

`kept_REAL` moves −3 (roster) / −10 (catalog) at max-removal (T=5,C=1). Attribution:

- **ROSTER −3 = 0 genuine paralogs lost.** All three are functionally-unrelated genes glued by an Alu that
  the cDNA-loose oracle itself mislabelled REAL: ASB6+NTMT1 (mult 503), EIF1AX+LOC (mult 497),
  ATP6V1C2+PDIA6 (mult 7). Cutting them is arguably correct.
- **CATALOG −10 contains only 2 genuine divergent cross-gene paralogs:** ZNF224 (KRAB-ZNF, exon-id 0.686)
  and RFPL3 (ret-finger, exon-id 0.639) — nicked by the strict 0.70 full-exon cutoff, **both outside the
  DNA-oracle denominator**. This is the irreducible EXON_CONTAINMENT 0.70-boundary limit, NOT a repeat-
  conjunct fault; the other 8 catalog losses are Alu/repeat-glued mislabels. Raising T shrinks the bleed.

Crucially, the conjunction does **not** bleed the way plain exon-containment did: the pure single-gene
fragmented reals that containment destroyed 1:1 are all spared (see below).

## What actually does the work (causal-role honesty)

Plain **(no-full-shared-exon AND same-gene guard, no repeat conjunct)** already cuts **32/35**. The repeat
conjunct **narrows** this to 24/35 (T5C1) / 15/35 (T8C2) by **abstaining** on the low-mult-glue tail. The
recall-safety on **fragmented single-gene reals** comes from the **same-gene GUARD, not the repeat
conjunct**: control fids **235** (max-cross-mult 20, m8=3) and **477** (mult-497 Alu, m8=8) both carry
strong repeat bridges yet are cut 0 because one gene spans the components (`guard_ok=0`). So:
**guard + no-full-shared-exon = the recall-safety; the repeat conjunct = a precision knob** (abstain vs cut
on the low-mult tail). The "multi-repeat-bridge" name credits the repeat conjunct more than its causal
share — stated plainly here.

---

## Connection to the shipped repeat-hub gate

`family_rna_refine.py` (DEFAULT-ON, `--no-repeat-gate` opt-out) cuts a between-gene edge iff a **single**
edge has `min_shared_mult ≥ 20` (`REPEAT_MULT_MIN = 20`, the extreme tail: 92.7 % RepeatMasker-concordant,
separating the fam17 hub from GSTM2 max-9 / MAGE max-8, both with zero edges ≥ 20). That is a
**single-extreme-EDGE** gate. The DISCONNECTED FPs **survive** it precisely because their bridges are
**sub-20 / distributed family-level**, not one extreme edge.

This gate is the **family-level multi-repeat extension** of that shipped gate: instead of one edge at
mult≥20 it counts **several** cross-component nodes at a **lower** threshold `T`, and — critically — makes
that safe by ANDing the **no-full-shared-exon** conjunct + **same-gene guard**. That conjunction is what
lets the multiplicity threshold drop from 20 to 5 **without** touching GSTM2/MAGE — the exact thing the
shipped single-edge gate could not do.

---

## Verdict

**YES-BUT (precision-safe, narrow, purity-positive).** The multi-repeat-bridge conjunction cleanly removes
the **repeat-bridged half** of the DISCONNECTED over-merge FP class (24/35 at T5C1, 15/35 at the C≥2
"several" point), improves E_p purity **+0.034**, holds the named R_oracle at **50/57** at every threshold,
and spares GSTM2/MAGE + all 23 pure single-gene controls at every point in both scopes. It does **NOT**
bleed into real single-gene paralogs the way plain exon-containment did — recall-safety comes from the
**conjunction + same-gene guard, not the multiplicity threshold** (so T can drop to 5). The other ~half of
DISCONNECTED (max-cross-mult ≤ 8) is **correctly abstained** — glued by low-mult short matches, not
high-copy repeats. The only residual is, catalog-wide, **2 cross-gene divergent paralogs** (ZNF224, RFPL3)
nicked at the strict 0.70 full-exon boundary — the irreducible exon-containment limit, outside the DNA
gold, shrinkable by raising T.

**Return line:** DISCONNECTED removed **24/35** (T5C1) / **15/35** (T8C2, the C≥2 "several" point) —
named above; GSTM2/MAGE **spared 0/5** (structurally, at every point, both scopes); **R_oracle 50/57
HELD** (set-identical, though insensitive-by-construction — `kept_REAL` is the operative probe: roster
−3 = 0 genuine paralogs, catalog −10 = 2 genuine near-0.70 paralogs); **E_p Δ = +0.034** (0.8694 → 0.9034).

## Recommendation: WIRE IT as a 3rd VG-native gate, DEFAULT-ON / opt-out

Wire the conjunction predicate as the **third VG-native family-refinement gate**, alongside the shipped
single-edge repeat-hub gate, **default-on with an opt-out** (mirroring `--no-repeat-gate`). Justification:
it removes the **dominant residual over-merge class**, improves purity, holds every recall-safety invariant
the audits asked about, and is VG-native (canonical-minimizer multiplicity + genome-base exon homology,
library-free). Two deployment forms:

1. **Recommended default (roster / conservative-catalog):** cut the DISCONNECTED-FP conjunction at
   **T=8, C=2** (the genuine "several-repeat" point) — 15/35 removed, P_Ep +0.022, kept_REAL −1, **0
   genuine paralogs lost**, R_oracle held. Zero-genuine-loss and honestly "several repeats".
2. **Aggressive (T=5, C=1):** 24/35 removed, P_Ep +0.034, still 0 genuine paralogs lost in roster scope.
   Use when maximum over-merge removal is wanted and the 2 near-0.70 catalog paralogs are acceptable
   (they fall outside the DNA gold and are the exon-containment 0.70-cutoff limit, not a gate fault).

Keep the **0.70 exon-id cutoff** as the one tunable to watch: relaxing it toward ~0.62 with the repeat
conjunct retained would recover ZNF224/RFPL3, at some precision cost — a follow-up, not a blocker. The gate
should NOT be shipped as multiplicity-alone (would shred GSTM2) nor exon-containment-alone (would shred
fragmented reals 1:1) — only the **conjunction + same-gene guard** is safe.
