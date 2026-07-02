# VG repeat catalog — a library-free, read-derived repeat definition in the variation-graph framework

**Code:** `bench/vg_repeat_catalog.py` · **Outputs:** `bench/vg_repeat_catalog.tsv` (215,302 node rows mult≥2 + 10,755 eval-edge rows), `bench/vg_repeat_catalog.json`
**Determinism:** re-run byte-identical — tsv `ed8ce6d906aeca3b1efd2da23df34761`, json `a73620dae82fca0407a0114625150609` (2-bit canonical encoding, leftmost-tie minimizers, `PYTHONHASHSEED=0`, sorted component folds).

## 1. The object — a canonical-minimizer MULTIPLICITY graph (library-free)

This is a novel O1 object parallel to the family definition: a repeat is defined **from graph structure over the de-novo transcripts alone**, with no RepBase / Dfam / RepeatMasker in the definition (soft-mask is used ONLY as external validation in §4).

- **Node** = a canonical (strand-agnostic) `(k=15, w=10)` minimizer value over the soft-masked spliced exonic sequence of each locus. Canonical (`min(fwd, revcomp)` 2-bit code) ⇒ an inverted paralog shares the same node; case-insensitive ⇒ RepeatMasker lowercase changes no node (library-free confirmed).
- **Universe** = the 3,462 gene-assigned loci that form E_r families (2,702 distinct genes).
- **Gene node-SET** = union of minimizer values over the gene's loci (order discarded → a bag of minimizers).
- **Path-multiplicity(node)** = # DISTINCT GENES whose node-set contains it. This is *set membership*, not k-mer occurrence count — the invariant `multiplicity ≤ n_transcripts` holds by construction (top node mult 503 < 588 transcripts).
- **Cycle flag(node)** = the minimizer recurs (≥2 runs) within one transcript (tandem/cyclic). This is the *only* place ordered-path structure is consulted.
- **REPEAT node** := `multiplicity ≥ m` (m swept) OR cyclic.
- **Repeat-aware edge oracle** := keep a family edge (A,B) iff the node-set intersection `A ∩ B` contains a node of multiplicity `< m` (a unique/low-mult minimizer); i.e. MASK repeat nodes before edge formation. Score(A,B) = min multiplicity over shared nodes (higher ⇒ over-merge).

**Graph:** 1,726,242 nodes; 215,302 shared (mult≥2); 6,310 cyclic.

> **Honest naming.** The delivered object is the *documented alternative* of the two the brief allowed — a canonical-minimizer / bag-of-minimizers **multiplicity** graph, not the maximal-shared-aligned-*segment* graph. The edge oracle is pure node-set intersection + min multiplicity; ordered "path" structure is used only for the cyclic flag. "Variation/segment graph" and "path-multiplicity" are therefore read as *minimizer-multiplicity / set-membership*, and this is stated in the code docstring.

## 2. Is VG-node-multiplicity a cleaner real-vs-overmerge separator than the scalar read-degree / soft-mask?

AUC(genuine-FP > TP) — higher = cleaner over-merge separator. **The only apples-to-apples comparison is on identical rows.**

| separator | AUC | rows (gen / TP) | note |
|---|---|---|---|
| **VG_min_shared_mult — matched (identical rows)** | **0.686** | 1773 / 416 | the honest number |
| core_recip — same rows | **0.814** | 1773 / 416 | beats VG |
| softmask — same rows | **0.774** | 1773 / 416 | beats VG |
| aln_frac — matched (sparse) | 0.853 | 1763 / 90 | beats VG |
| VG_min_shared_mult — full VG-covered pop | 0.774 | 3727 / 449 | **inflated, see below** |
| VG out-of-fold (5-fold CV, 670 comps, seed 0) | 0.742 | — | mean test AUC |
| prior read-multimap degree | 0.64–0.69 | — | the "scalar read-degree" |
| prior homology-graph degree | 0.83 | — | not beaten |

**Honest read (corrects the earlier favorable spin):**
- On **identical-methodology rows** VG = **0.686**, which is **worse than soft-mask (0.774)**, worse than core_recip (0.814), and worse than aln_frac (0.853). It does **not** "tie soft-mask."
- Versus the scalar read-degree, VG only **ties the top of the read-multimap range (0.686 ≈ 0.69)** on identical rows — it does not cleanly beat it, and it does not beat the homology-graph degree (0.83).
- The apparent VG-full **0.774** is **inflated**: the extra rows it includes over the matched set are **1,954 genuine / 33 TP (98.3% genuine)**, a near-degenerate easy subpopulation on which VG scores only 0.725. VG-full 0.774 vs soft-mask-matched 0.774 is apples-to-oranges, not a tie.
- Out-of-fold CV = 0.742, still below core (0.814), aln_frac (0.853), homology-degree (0.83). VG does **not** move the P/R frontier.

**Conclusion for §2:** as a *general scalar over-merge separator*, VG-node-multiplicity is **the same soft signal in graph form** — cleaner than the raw read-degree but not cleaner than the alignment signals.

## 3. fam17 hub cut, GSTM2 / MAGE spared, and the recall cost (m = 5 operating point)

Where the VG object *does* earn its keep is as a **targeted, library-free repeat identifier**, not as a global separator.

| case | genes | edges | cut@m5 | frac cut | median min-shared-mult | shared-node max mult | hi-mult softmask |
|---|---|---|---|---|---|---|---|
| **fam17** (16-family repeat-bridge hub) | 27 | 301 | **252** | **0.837** | 11 (max 38) | **503**, cyclic | **0.999** |
| GSTM2 (protein DOMAIN — neg. control) | 20 | 78 | 28 | 0.359 (64% kept) | 4 (max 9) | 503* | 0.315 |
| MAGE (cardinality — neg. control) | 11 | 55 | 16 | 0.291 (71% kept) | 4 (max 8) | 10 | 0.112 |

- **fam17: clean WIN.** The bridge is a genuine library-free repeat node — multiplicity up to **503 distinct genes**, cyclic, and its high-mult shared nodes are **99.86% soft-masked**. The oracle shatters **83.7%** of the 16-family hub.
- **GSTM2 and MAGE: majority-SPARED (64% / 71% kept).** Their edges survive because the *min* shared multiplicity per edge is low (median 4), i.e. each pair also shares a unique coding node; their high-mult shared nodes are mostly unmasked coding domain (softmask 0.31 / 0.11), correctly **not** flagged as repeats. (GSTM2's shared-node *set* does touch mult-503 Alu nodes, but those are masked away and the surviving edges rest on low-mult coding nodes.)
- The discriminator that cleanly separates the true hub from the two controls is **multiplicity magnitude (503 vs 4–10) × softmask (0.999 vs 0.31 / 0.11)** — the true repeat is both.

**Recall cost (honest — the boundary is soft):** at m=5 the pure-multiplicity threshold is not free. It keeps **88.7% of TP edges** and cuts **59.5% of genuine over-merges**, but it also cuts **≈54% of truth-bar borderline-real paralogs** (truthbar recall 0.461) and loses **36% / 29%** of GSTM2 / MAGE edges as collateral. Out-of-fold CV: TP-recall 0.849, genuine-cut 0.566. **It does cut real families.** A softmask-*gated* repeat definition (fam17 0.999 vs GSTM2 0.31) would remove almost all of this collateral — but softmask is validation, not library-free.

## 4. RepeatMasker concordance (validation only — softmask absent from the definition)

Do high-path-multiplicity nodes correspond to actual soft-masked genome?

- Node-level AUC(mult → soft-masked ≥ 0.5) = **0.574** overall (weak; mult-2 nodes are mostly unique coding, baseline 20.6% masked). Spearman(mult, softmask) = 0.179.
- **Strongly monotone at the tail:** mult≥2 → 23.8% masked · mult≥3 → 35.3% · mult≥5 → 48.5% · mult≥10 → 74.4% · **mult≥20 → 92.7%**. The top nodes (mult 503, 497, 452 …) are all softmask = 1.0.

**Honest caveat.** The extreme high-mult tail is a **mix**: a few of the very top nodes are **low-complexity** poly-A minimizers (2-bit codes 0/2/8 = `AAAAAAAAAAAAAAA…`), trivially both high-mult and masked; but most are **genuine interspersed-TE fragments** — the top nodes `AAAGTGCTGGGATTA`, `AATCCCAGCACTTTG`, `AGGCTGAGGCAGGAG` are canonical **AluY consensus** motifs. So the 93% high-tail concordance reflects **both** low-complexity **and** real interspersed TE rediscovered library-free — not purely one or the other. (RepeatMasker soft-mask itself includes low-complexity, so this is fair for a validation.)

## 5. Honest verdict

The VG repeat catalog is a **worthwhile, defensible VG contribution as a library-free repeat *object* — but not as a new precision separator that beats alignment.** Framed properly (nodes = canonical minimizers, gene = a set of nodes, multiplicity = distinct-gene set-membership, cycle = intra-transcript tandem), high-multiplicity/cyclic nodes constitute a genuinely library-free repeat catalog that **rediscovers RepeatMasker repeats at the high-multiplicity tail (mult≥20 → 93% soft-masked, top nodes = AluY consensus + poly-A) with no external library**, and that **surgically identifies the true fam17 repeat-bridge hub** (mult-503, 99.9%-masked, cyclic bridge node; 84% of the 16-family hub cut) while **majority-sparing** the GSTM2 protein domain (64% kept, softmask 0.31) and the MAGE cardinality FP (71% kept, softmask 0.11) — both honestly low-multiplicity, low-softmask, and therefore *not* flagged as repeats. As a general real-vs-overmerge separator, however, it is **the same soft signal in graph form**: on identical rows VG multiplicity (AUC **0.686**) is *worse* than soft-mask (0.774), core_recip (0.814), and aln_frac (0.853), only *ties the top of the read-multimap-degree range* (≈0.69), does not beat homology-graph degree (0.83), does not move the P/R frontier, and at the m=5 operating point cuts ≈54% of borderline-real paralogs — so pure multiplicity-masking is not a free precision lever. The genuinely clean lever remains **multiplicity gated by coding/softmask** — the same conclusion as the prior TE work, now expressed as a defensible variation-graph object (nodes / paths / multiplicity / cycles) rather than a separator that beats alignment.

**Headline numbers.** VG graph 1,726,242 nodes / 215,302 shared / 6,310 cyclic over 3,462 loci / 2,702 genes · separator AUC(genuine>TP) identical-rows **0.686** vs soft-mask 0.774 / core_recip 0.814 / aln_frac 0.853 / read-degree 0.64–0.69 / homology-degree 0.83 (VG-full 0.774 inflated by 1,954-gen/33-TP rows) · CV out-of-fold AUC 0.742, TP-recall 0.849, genuine-cut 0.566, truthbar-recall 0.461 · **fam17 252/301 edges cut (83.7%), bridge mult 503 cyclic 99.9% masked** · GSTM2 64% kept (softmask 0.31) · MAGE 71% kept (softmask 0.11) · RepeatMasker concordance mult≥20 → 92.7% masked (AUC 0.574, spearman 0.179).
