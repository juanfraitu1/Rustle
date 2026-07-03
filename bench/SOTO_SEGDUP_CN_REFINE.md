# Soto-style DNA-calibrated refine: does the SEDEF-segdup(+famCN) DNA signal split the domain-share FP residual from real families?

**DNA-INFORMED, NOT the RNA-only gate.** This is the RNA-catalog analog of Soto et al. 2025's family-cleaning:
Soto define a family as **SD98** (segmental duplications at >=98% identity) + keep exon matches covering >99% of
each exon + **famCN** refine (group only genes whose read-depth copy-number mean-abs-deviation < 1). Both the SD98
leg and the famCN leg are **DNA**. We stand in the gorilla analogs — the **SEDEF segdup map**
(`/mnt/c/Users/jfris/Desktop/final.bed`, 253,029 pairs, SD98 analog) and the **mGorGor1 diploid CN**
(`bench/diploid_cn_oracle.tsv`, 51-gene multi-copy famCN analog) — and ask whether the DNA duplication signal
cleanly separates the **domain-share FP residual** (functionally-unrelated genes sharing one exon: RHD+SDHD,
RBBP4+SYNC, DEDD+NIT1, DIMT1+KIF2A, …) from real paralog families, and **at what recall cost on DIVERGENT real
paralogs**. Deterministic (PYTHONHASHSEED=0; re-run JSON+TSV diff=0, md5 json `d37a1dcd`, tsv `92fdfcb4`).

This is labeled a **Soto-style DNA-calibrated REFINE / validation, like the diploid oracle** — it is NOT the
RNA-only edge rule. Substrate = the same 5,571 direct multi-exon E_r edges as the RNA-only colinear gate
(416 TP + 3,382 truthbar + 1,773 genuine), so segdup-refine and colinear-gate are apples-to-apples.

Artifacts: `bench/soto_segdup_cn_refine.py` / `.json` / `.tsv`. SEDEF overlap uses the actual segdup **PAIR**
(the two loci sit on the two ends of ONE final.bed line = locus A duplicated to locus B), not mere co-location;
identity = fracMatch (col 21, 1-indexed), fracMatchIndel (col 22) carried. Loci -> coords via
`denovo_transcripts.meta.tsv`, 0 coord drops.

---

## (1) DOMAIN-SHARE FP removed vs real kept — PARTIAL separation, with two honest failures

### DNA-reachable domain-shares (shared exon is NOT itself a duplication) — cleanly removed
SD90 cleanly cuts the **5 named single-exon domain-shares** whose shared exon is not a segdup — all at
**sd = 0.000** (no cross-pair in SEDEF at all, or both loci in segdups but no A->B pair):

| pair | cross_n | sd (cross fracMatch) | why removed |
|---|---|---|---|
| MOV10+RHOC | 0 | 0.000 | **no segdup at all** (pure case: a_in=b_in=0) |
| RBBP4+SYNC | 0 | 0.000 | both loci in segdups but not to each other |
| DEDD+NIT1 | 0 | 0.000 | no A->B pair |
| ABITRAM+CTNNAL1 | 0 | 0.000 | both loci self-align to a THIRD locus, same end -> no cross-pair |
| DIMT1+KIF2A | 0 | 0.000 | no A->B pair |

At the pair level the true single-exon target (**ep_impure genuine, colinear==1, n=15**) is only **20%
cross-supported at SD90 -> 80% removed** (12/15). This is the DNA-reachable core of the FP residual.

### DNA-UNreachable domain-shares (shared block IS a segdup) — survive SD98
Two named "domain-shares" resist because the shared block is itself a genuine duplication — SD98 alone cannot
tell them from a real recent paralog:

| pair | cross_n | sd | colinear | verified in raw final.bed | SD-floor that removes it |
|---|---|---|---|---|---|
| **RHD+SDHD** | 2 | 0.977 | 1 | cross-chrom NC_073224.2 <-> NC_073233.2, ~1 kb, fm 0.964/0.977 | only **SD98** |
| **LOC129532935+SORT1** | 1 | 0.987 | 8 | 77 kb cross-pair NC_073224.2:116.67M <-> 130.66M, fm 0.9866 | **none** (survives SD98) |

### The full ep_impure roster is impure — worse than divergent reals
The full **ep_impure genuine roster (n=58)** is **58.6% cross-supported at >=0.90** — WORSE separation than divergent
reals — because the roster is impure: its multi-exon members are genuine high-identity segdups. Pair-level removal:
**41% (24/58) at SD90, 59% (34/58) at SD98.** So the clean col==1 subset removes at 80%, but the mixed roster only
41–59%.

### Reals are cleanly kept in the near-identical regime
REAL **near-identical (TP, n=416)** sit cleanly IN the SEDEF signal: cross-support **>=0.90 = 84.4%**
(any-id 92.3%, both-in-segdup 94.7%, >=0.98 50.7%). The SD98 leg keeps recent paralogs.

---

## (2) DECISIVE — DIVERGENT-real recall cost: fraction of divergent paralogs captured in SEDEF

**truthbar (n=3,382)** = E_p-protein-confirmed DIVERGENT real paralogs that DNA did not confirm. Fraction whose
loci are a SEDEF segdup PAIR at the SD floor:

| SEDEF cross-pair floor | captured | CUT (recall cost) |
|---|---|---|
| **>=0.90 (SD90)** | **2.10% (71/3382)** | **97.90%** |
| >=0.95 | 0.71% (24) | 99.29% |
| >=0.98 (SD98) | 0.33% (11) | 99.67% |
| any-identity cross-pair | 18.63% (630) | 81.37% |
| each-locus-in-any-segdup (weak) | 66.50% | — |

**The SD90/SD98 restriction CUTS ~98% of divergent reals.** This is not a threshold artifact — dropping the floor
0.98 -> 0.90 recovers only 0.33% -> 2.10%; it is **erosion**, not a knob.

**Two distinct erosion drivers** (both disclosed in the JSON; splitting them as the reviewers recommend):
- **No SEDEF pair at all: 81.4% (2,752/3,382)** — the divergent paralog does not self-align at any identity.
- **Pair present but < 0.90: 16.5% (559/3,382)** — a cross-pair exists but the shared block has diverged below
  the SD floor. (Of the 630 truthbar with any cross-pair, 559 = 88.7% are < 0.90; only 11 reach >=0.98.)

**Monotone erosion along the RNA divergence axis** — the DNA erosion confound mirrors the RNA one exactly.
Fraction captured at SEDEF >=0.90, stratified by mean-exon RNA identity:

| RNA mean-exon id | n | captured at SEDEF >=0.90 |
|---|---|---|
| >=0.95 | 29 | 65.52% |
| 0.90–0.95 | 81 | 19.75% |
| 0.80–0.90 | 1,207 | 2.82% |
| 0.70–0.80 | 1,423 | 0.00% |
| < 0.70 | 642 | 0.31% |

The ~98% cut, if anything, **overstates** cost (truthbar is assumed 100% real) and never understates it.

---

## Soto-style REFINE applied to the shipped 573-family catalog

`refine_stage()` mirrors the RNA-only `colinear_multiexon_gate.gate_stage` machinery EXACTLY (locus subgraph per
catalog family + same-gene backbone guard -> re-derive components -> drop <2-loci -> re-eval via shipped
`eval_partition`/`oracle_residuals`/gene-projection relabel), but the KEEP rule is Soto's: **keep a cross-gene
family edge iff the gene-pair is SEDEF-segdup-supported at fracMatch >= floor** (sweep SD90/SD95/SD98), plus an
isolated famCN (diploid-CN) leg. Re-derive baseline (no cut, matches the colinear gate):
**572 fam / E_p 0.8916 / R_oracle 50/57**, 1,812 co-membered substrate pairs.

| floor | fam | E_p | E_p delta | R_oracle | oracle lost | ep-impure FP split | REAL split (TP + truthbar) | split-FP precision | GSTM2 / MAGE |
|---|---|---|---|---|---|---|---|---|---|
| **SD90** | 409 | 0.912 | **+0.0204** | **49/57** | LOC101127437 | 20 | 735 (52 + 683; 638 at id<0.90) | **30.7%** | spared / spared |
| SD95 | 370 | 0.9108 | +0.0192 | 48/57 | +LOC115931067 | 25 | 817 (96 + 721) | 31.0% | SPLIT / spared |
| SD98 | 329 | 0.924 | **+0.0324** | **43/57** | +MAGE 129529978/986 (7 total) | 31 | 894 (163 + 731) | 32.0% | SPLIT / SPLIT |

**No SEDEF floor holds R_oracle >= baseline (50/57)** -> `best_floor = None`. Even the least-damaging floor (SD90)
loses one oracle gene and is **strictly dominated** by the RNA-only colinear K=0 gate (which holds 50/57 at
**72.7% split-FP precision** vs segdup SD90's **30.7%** — the segdup refine cuts >2 real per FP removed, a far
blunter instrument). Every E_p gain is bought by breaking R_oracle and cutting divergent reals.

Wording precision (honest): removing the **5 sd = 0.000** single-exon FP edges by itself costs **no** oracle genes,
but **SD90 as an operating point is 49/57** — that one-gene loss (LOC101127437) is DIVERGENT-real collateral, not
those 5 edges. "Removes the clean FP at zero oracle cost" is true only of the 5 edges, never of SD90 as a whole.

### famCN (diploid-CN) leg — too sparse to adjudicate cardinality
DNA-reachable for only **9/5,571 pairs** (both genes in the 51-gene multi-copy oracle). It flags **6 CN-inconsistent
gene-pairs** — counting precisely: **4 GSTM2-anchored** (CN19 vs 3/13/16/4 = LOC101129940/LOC115930164/
LOC115930576/LOC115935025) **+ 1 GSTM-cluster TP FALSE-FLAG** (LOC115930164+LOC115930576, CN13 vs 16 — MAD<1 is
too strict for a large amplicon) = **5 GSTM-family pairs**, **plus FOXO1+LOC115933254** (CN2 vs 8) = **6 total.**
famCN-only splits just 1 genuine FP (FOXO1) and drops R_oracle to 49/57 (same FOXO1 casualty as colinear K=1).
Critically, cutting the CN-inconsistent GSTM2 edges leaves **GSTM2 (5 fam) and MAGE (2 fam) cardinality UNCHANGED**
— the amplicons stay connected via the same-gene backbone. The famCN leg is too sparse to touch cardinality: the
domain-share partners are single-copy genes ABSENT from a multi-copy oracle (that absence is itself weakly
consistent with "not a duplicated family").

---

## Verdict — how much of the ~8–11% E_p FP residual is DNA-REACHABLE vs irreducible

- **DNA-REACHABLE** = clean single-exon domain-shares whose shared exon is **NOT** itself a duplication
  (MOV10+RHOC, RBBP4+SYNC, DEDD+NIT1, ABITRAM+CTNNAL1, DIMT1+KIF2A; ~20 ep-impure edges at sd = 0.000). SD90
  removes these — the ~80% of the col==1 single-exon target — and the removal of those specific edges costs no
  oracle genes.
- **DNA-UNREACHABLE / irreducible** = "shared-domain-IS-a-segdup" cases (**RHD+SDHD** sd 0.977, **LOC129532935+SORT1**
  sd 0.987, colinear=8). SD98 removes only RHD+SDHD, and only by dropping into the near-identical regime — which
  destroys R_oracle (43/57) and MAGE cardinality; LOC129532935+SORT1 survives even SD98. **SD98 alone cannot
  separate a domain-share whose shared block is itself a >=98% duplication from a real recent paralog.**

The DNA duplication signal is a **validator for the recent/near-identical regime only.** It cleanly keeps
near-identical reals (84%) and cleanly removes genuine single-exon domain-shares (80%), but at a **catastrophic and
near-irreducible recall cost on DIVERGENT reals: ~98% eroded below SEDEF's ~90% floor** (81% have no SEDEF pair at
all; 16.5% have a pair that has diverged below the floor). It imposes the **same divergence-erosion ceiling at the
DNA level that the RNA gate hits**, it is a **blunter instrument** than the RNA-only colinear K=0 gate (30.7% vs
72.7% split precision, 49/57 vs 50/57), and its famCN leg is **too sparse (9/5,571)** to touch cardinality.
**It cannot be the primary family-definition criterion for divergent paralogs.**

### Connection to Soto 2025
Soto needed **SD98 + >99%-exon-coverage + famCN — all DNA** — precisely because the domain-share is **DNA-bound**:
the RNA/exon signal alone cannot separate a functional domain-share from a paralog (both share a full exon), so
Soto reach for the DNA duplication signal. Our result **confirms the domain-share is DNA-bound** (the RNA-only gate
cannot fully separate it), but also shows the DNA signal is **not a complete fix**: (a) Soto's SD98 is engineered
for the **recent/near-identical** SD regime and imposes the same divergence-erosion ceiling we hit — 98% of
divergent gorilla reals fall below the SD90 pairwise floor; and (b) the "shared-domain-IS-a-segdup" residual
(RHD+SDHD) survives even SD98, which is exactly why Soto lean on the additional >99%-exon-coverage + famCN legs to
finish. Our famCN analog (mGorGor1 diploid CN) is present but too sparse to close that gap on this RNA substrate.

### Recommendation — VALIDATION-ONLY, not a wired edge-cutter
**Keep the Soto refine as a DNA-calibrated VALIDATION overlay (like the diploid oracle), do NOT wire it as a
family-definition edge-cutter — not even as an opt-in `--soto-refine`.** Rationale: no SEDEF floor is recall-safe
(best_floor = None), and as an operating point it is strictly dominated by the RNA-only colinear K=0 gate on the
identical task (better R_oracle, better split precision). If any exposure is wanted, expose the **per-edge
segdup-support as a read-only DNA-CONFIRMATION annotation** (a cross-gene edge whose gene-pair is a >=0.90 SEDEF
pair = "DNA-confirmed"; TP near-identical edges are 84% DNA-confirmed, divergent edges 2%), which tags the
near-identical/DNA-confirmed subset **without cutting any edge** — never as a default or opt-in cutter. This keeps
the honest boundary: the DNA signal validates the recent regime and flags the DNA-reachable single-exon
domain-shares, but the divergent-paralog family definition must remain RNA-homology-driven.

---

Deterministic re-run: `PYTHONHASHSEED=0 python bench/soto_segdup_cn_refine.py` -> JSON+TSV diff=0
(md5 json `d37a1dcd`, tsv `92fdfcb4`).
