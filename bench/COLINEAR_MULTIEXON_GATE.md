# Colinear multi-exon family gate — does "require ≥2 shared exons in colinear order" remove the DOMAIN-SHARE FP while sparing real paralogs?

**2026-07-03.** Tested the user hypothesis for killing the **DOMAIN-SHARE FP residual** (E_p-impure pairs that the shipped edge rule `core_recip ≥ 0.19 AND aln_frac ≥ 0.24` glues into families because `aln_frac` = the SINGLE longest shared spliced-exon block / shorter transcript — it requires neither ≥2 shared exons NOR the same exon order).

**User model of reality:** a recent paralog diverged from a *whole shared gene* and retains **≥2 exons in colinear order**; a domain-share is **one isolated exon** in different flanking context. A handful of exons may reorder (gene conversion / exon shuffle), but it is unlikely that *all but one* do. **Proposed gate:** within each family, cut a cross-gene edge iff `colinear_shared_exons ≤ K` (K=1 ⟺ "require ≥2 colinear exons"), same-gene guard, drop <2-loci.

Script `bench/colinear_multiexon_gate.py` (+ `.tsv` / `.json` / `colinear_multiexon_cache.tsv`). Deterministic, PYTHONHASHSEED=0, RNA-only (exon structure + genome-base exon homology; DNA / E_p = validation only). Byte-identical on rerun (md5 `e3f35e4b` / `bcc6f642` / `505793b2`). Cross-checked against the shipped `family_level_pr_current` machinery: gate baseline reproduces catalog **573 fam / E_p 0.8918 / R_oracle 50/57** exactly.

## Method

**`colinear_shared_exons`** for a gene-pair (dn_a, dn_b): reconstruct per-exon spliced sequences (`denovo_skeletons` introns → exon complement, transcript 5'→3' per strand — identical machinery to `recombination_bridge_detector.family_exons` / `exon_match_tensor`), align every exon of one against every exon of the other (edlib HW infix, fwd+revcomp, **id ≥ 0.70**), take each exon's best target, then **LIS = the longest strictly order-preserving chain of DISTINCT matched target exons**. Max over the two directions (recall-safe). It is a raw **count**, not a fraction: an isolated shared domain → 1; a shared multi-exon gene body → ≥2.

**Substrate** = the 5571 MULTI-exon core-present direct E_r edges in `bench/exon_containment_criterion.tsv` (the shipped oracle's learn population, representative loci already selected): **416 TP** (DNA-confirmed, near-identical), **3382 truthbar** (E_p-protein-confirmed **DIVERGENT** real paralogs = the confound), **1773 genuine** (over-merge FP = the domain-share residual). REAL paralogs = TP + truthbar; the divergent reals = truthbar.

**Gate stage** operationalizes faithfully to shipped `family_level_pr_current`: per catalog family, locus subgraph of the shipped E_r edge graph + a same-gene backbone chain (the guard — same-gene loci never separate); cut every cross-gene edge with gene-pair `colinear ≤ K`; re-derive components; drop <2-loci; re-eval via shipped `eval_partition` + `oracle_residuals` + gene-projection relabel. Deltas are measured against the SAME operationalization with no cut (re-derive baseline **572 / E_p 0.8916 / R_oracle 50/57**; the 1-family / 2-locus gap vs the 573 catalog is a re-derivation artifact, NOT the gate). Colinear covers 1812/1817 intra-family cross-gene gene-pairs.

---

## ⭐ DECISIVE evidence — divergent-real-paralog colinear-exon distribution (at the specified id ≥ 0.70)

| class | n | col=0 | col=1 | **col≥2** | col≥3 | median |
|---|---|---|---|---|---|---|
| TP (near-identical real) | 416 | 7.2% | 9.4% | **83.4%** | 74.3% | 4 |
| **truthbar (DIVERGENT real)** | 3382 | 19.0% | 34.7% | **46.3%** | 23.5% | **1** |
| genuine (domain-share FP) | 1773 | 45.5% | 20.1% | 34.4% | 25.6% | 1 |

**Only 46.3% of divergent (truthbar) real paralogs retain ≥2 colinear exons at id ≥ 0.70.** 34.7% erode to a single domain (col=1) and 19.0% share no exon ≥0.70 at all (col=0). **Median colinear = 1.** So the majority (53.7%) of divergent reals are **CUT at K=1** — the user's structural claim does NOT hold at the specified floor.

Stratified by mean matched-exon identity: id≥0.95 → 79% keep ≥2; 0.90–0.95 → 74%; 0.80–0.90 → 58%; 0.70–0.80 → 55%; "no exon ≥0.70" (n=642) → 0%. Real paralogs with mean-exon-id <0.90: **45.3% ≥2** (n=3342); <0.80: **37.2% ≥2** (n=2099).

### Root cause — nucleotide-identity-FLOOR artifact, not lost structure

col≥2 % per class as the exon-match identity floor varies:

| floor | TP | **truthbar (divergent real)** | genuine FP | real–FP gap |
|---|---|---|---|---|
| id≥0.50 | 93.8 | **95.4** | 82.2 | 13 pt |
| id≥0.60 | 84.6 | **74.6** | 45.8 | **29 pt** |
| id≥0.65 | 83.9 | 57.9 | 37.6 | 20 pt |
| **id≥0.70** | 83.4 | **46.3** | 34.4 | 12 pt |
| id≥0.80 | 82.2 | 21.6 | 27.5 | −6 pt |

**The colinear STRUCTURE is real in divergent paralogs — 95.4% keep ≥2 at id≥0.50** — but their exons sit at ~0.55–0.70 nt identity and fall below the 0.70 floor, so they *present* as isolated domains, indistinguishable from single-domain FPs. Verified end-to-end: of the 1815 divergent reals cut at K=1, **1660 (91.5%) recover col≥2 once the exon-match floor drops to id≥0.50**. Any nucleotide-identity-thresholded exon criterion (containment or colinear-count) inherits this. Class separation is *weakest exactly at 0.70* (truthbar 46% vs genuine 34% = 12 pt); the best operating point is ~id≥0.60 (divergent-real recall 74.6% vs FP 45.8%, **29 pt gap**) — better, but even there 25% of divergent reals are cut and 46% of FPs kept.

---

## Gate outcome on the shipped catalog (best K, R_oracle, E_p, splits)

| K (gate) | n_fam | E_p (Δ vs baseline) | R_oracle | FP split (ep-impure) | REAL split (TP + divergent) | split-FP-frac |
|---|---|---|---|---|---|---|
| **K=0** (cut col=0; zero shared exon) | 499 | 0.9078 (**+0.016**) | **50/57 HELD** | 72 (11) | 27 = 16 + 11 | 72.7% |
| **K=1** (require ≥2 colinear = USER GATE) | 428 | 0.9346 (**+0.043**) | **49/57 FAIL** (loses FOXO1) | 164 (26) | 67 = 40 + 27 (24 id<0.90) | 71.0% |
| K=2 (require ≥3) | 404 | 0.9332 (+0.041) | 46/57 FAIL (FOXO1 +3) | 222 (30) | 162 = 72 + 90 | 57.8% |

**Best K = 0** — the ONLY K that holds R_oracle 50/57 (recall-safe), E_p +0.016. But K=0 cuts only zero-shared-exon edges; it does **NOT** implement the user's "require ≥2 colinear" proposal (that is K=1). Per the requested breakdown:

- **(a) DOMAIN-SHARE FP removed.** K=1 splits all 5 present named exemplars — **MOV10+RHOC (col=0), RBBP4+SYNC (col=1), RHD+SDHD (col=1, aln_frac=0.91), DEDD+NIT1 (col=1), ABITRAM+CTNNAL1 (col=1, aln_frac=0.76)** — plus 21 more ep-impure pairs (26 total), and correctly SPARES **LOC129532935+SORT1** (col=8, structurally multi-exon, likely mislabeled FP). These are exactly the FPs the shipped `aln_frac` rule keeps. (KMT2C+TMEM128 has no direct multi-exon E_r edge in the substrate — honestly absent.) At K=0 only MOV10+RHOC (col=0) is removed; the conserved col=1 domain-shares SURVIVE.
- **(b) REAL paralog cost (divergent broken out).** K=1 splits **67 real pairs: 40 TP + 27 divergent-truthbar** (24 at mean-exon-id <0.90). Named divergent casualties are the nucleotide-floor col≤1 families sharing no exon at id≥0.70: **the whole KRAB-ZNF cluster** (ZNF41/ZNF585A/ZNF585B/ZNF84 …), **protocadherin-beta** (PCDHB2/5/15/16), **HNRNPH1+HNRNPH2, RPL10+RPL10L, HSPA1L+HSPA2, DCAF4+DCAF4L1, IFITM1+IFITM2, APOL2/3/4, NR2F1+NR2F2** — real paralog families, split. At the per-edge level K=1 cuts 1163 FP but 1815 divergent reals + 69 TP: **only 38.2% of cut edges are FP, 61.8% are real**; the col=1 bucket is **3.4:1 REAL:FP** (1212 real : 356 FP). Family-level co-membership (same-gene backbone + alternative paths) softens this to 71% split-FP at K=1, but does not rescue it.
- **(c) R_oracle 50/57.** HELD only at K=0. FAILS at K=1 (49/57, loses **FOXO1**) and K=2 (46/57). The FOXO1 casualty is a single catalog locus in the ANKRD18 segdup family; the gate cuts its ≤1-colinear cross-gene edges, isolating it to a singleton the <2-loci floor drops. (Honest nuance: FOXO1+LOC115933254 is itself a named multifam over-merge, so the metric is partly penalizing a correct FP removal — but per the shipped oracle, recall does not hold.)
- **(d) E_p purity delta.** +0.016 (K=0), **+0.043 (K=1)**, +0.041 (K=2). K=2 buys no purity over K=1 while doubling the real cost.
- **(e) GSTM2 / MAGE spared.** YES at all K. GSTM2 median colinear 8, 77% ≥2, only 23% cut (in 5 families, genes unchanged); MAGE median 2, 88% ≥2, 12% cut (UTR exons save the intronless retrocopies, 2 families). Cardinality preserved.

---

## Connection to the exon-containment falsification (`EXON_CONTAINMENT_CRITERION.md`)

The colinear-ORDER refinement was proposed precisely to succeed where **coverage-containment failed** (containment cut real paralogs ~1:1 because divergent paralogs share few exons; AUC 0.900 < aln_frac 0.915; could not cut the GSTM2 domain-hub). **It does not.** Adding order to the exon-match count changes *which* exemplars are cut (it now cleanly removes the *conserved* single-domain-share exemplars RBBP4+SYNC / RHD+SDHD / DEDD+NIT1 / ABITRAM+CTNNAL1 that containment left in, and it correctly keeps GSTM2/MAGE intact via their high colinear count), **but it inherits the identical failure mode on divergent paralogs** — because both criteria are gated on the same id≥0.70 nucleotide floor, and divergent real exons sit at 0.55–0.70. The order constraint is not the bottleneck; the identity floor is. This is exon-containment redux at the recall-critical divergent tail.

---

## Verdict

**The user's hypothesis FAILS at the specified id≥0.70.** Requiring ≥2 shared exons in colinear order (K=1) cleanly removes the *conserved* single-domain-share FPs (MOV10+RHOC, RBBP4+SYNC, RHD+SDHD, DEDD+NIT1, ABITRAM+CTNNAL1 — all 5), spares GSTM2/MAGE and cardinality, and lifts E_p by **+0.043** — but it **breaks R_oracle (50→49, loses FOXO1)** and **splits 67 real paralog pairs, 27 of them divergent** (the whole KRAB-ZNF / PCDHB / HNRNPH / RPL10L / HSPA / APOL / NR2F clusters). Only **46.3%** of divergent real paralogs retain ≥2 colinear exons at id≥0.70; median = 1; the confound wins roughly 1:1-to-worse with FP, exactly like plain exon-containment. The ≥2-colinear *structure* genuinely exists (95.4% keep ≥2 at id≥0.50, 91.5% of the cut divergent reals recover it), but the 0.70 nucleotide floor hides it, so divergent paralogs read as isolated domains.

### Recommendation — KEEP OUT (do not ship the user's gate as a 4th default-on gate)

- **Do NOT wire K=1 (≥2-colinear at id≥0.70) as a family gate.** It fails the R_oracle recall bar and cuts real divergent paralog families the pipeline claims elsewhere to recover.
- **The only R_oracle-safe operating point is K=0** (cut only ZERO-shared-colinear-exon edges: R_oracle held 50/57, E_p +0.016, 72.7% of splits are FP, removes MOV10+RHOC). If any colinear gate ships, ship K=0 default-on/opt-out — but note it leaves the col=1 conserved domain-shares (RBBP4+SYNC, RHD+SDHD, DEDD+NIT1, ABITRAM+CTNNAL1) in place and therefore does **not** deliver the user's hypothesis.
- **Recall-safe use of the user's ≥2 idea would require lowering the exon-match floor to ~id≥0.60** (divergent-real recall recovers to ~75% at a 29 pt FP gap) rather than 0.70 — better, but still not clean (25% of divergent reals cut, 46% of FPs kept).
- The residual domain-share FP is fundamentally a **DNA-family-boundary / cardinality** problem (finer divergence than RNA exon identity resolves), consistent with the containment finding — it is **DNA-bound**, not an RNA exon-structure gate away.
