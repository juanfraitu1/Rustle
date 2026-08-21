# O1 — census of incorrectly-called families: does the definition survive?

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

**Status 2026-08-15. Answer: the definition does NOT come through clean.** There are ~30 definitional
failures, and they are **not 30 separate problems — they are one mechanism with 30 instances.**

⚠ **INCOMPLETE RUN.** A session limit killed 9 of 25 agents: 5 classification chunks, all 3 adversarial
attacks, and the report agent. **138 cases were harvested; 105 were classified; 33 are unclassified.**
The three attacks — including the one aimed at whether the definitional count is *under*-reported —
never ran. **Treat 30 as a floor, not a final count.** Resume:
`Workflow({scriptPath: ".../o1-error-case-census-wf_49800d22-fa1.js", resumeFromRunId: "wf_49800d22-fa1"})`

---

## 1. The answer in numbers

| | count | note |
|---|---|---|
| cases harvested | **138** | memory + docs + empirical sweep of the shipped catalogs |
| classified | **105** | 33 unclassified — session limit |
| **DEFINITIONAL failures** | **30** | the rule admitted what it should not have. **These threaten O1.** |
| node-construction failures | 47 | the rule applied to a bad node. Definition innocent; catalog still defective. |
| not an error | 26 | correctly called after checking |
| unresolvable | 2 | sequence not recoverable |

**Measured rates, species kept separate, unit stated:**
* false-MERGE, **human**, unit = window: **2/150 = 1.33% [0.37, 4.73]**, power measured.
* false-MERGE against external truth (HGNC gene_group), **human**, unit = gene pair:
  **1/15 = 6.67% [1.19, 29.82]**.
* false-OMISSION, **gorilla**, unit = family: **9/162 = 5.6% [3.0, 10.2]**, upper bound 24.1%.

## 2. ⭐⭐⭐ The definitional failures are ONE mechanism: the min-length coverage denominator

`E_r` requires coverage ≥ 0.50 **of the shorter sequence**. That makes coverage **scale-free**: a
~1 kb dispersed repeat is ≥ 0.50 of *any* node under ~2 kb. Identity never binds — the offending edges
sit at **0.749–0.803**, far above the 0.60 floor. **It is the coverage clause, on a short node, that
admits them.**

**Exposed population is large: 352/1415 = 24.88% of shipped gorilla copies have reps ≤ 2,000 bp.**
⭐ **Re-measured 2026-08-20 on the current 627-family catalog: 432/2019 = 21.40%.** The mechanism
claim survives — still about one copy in five — but **quote 0.2140, not 0.2488**.
Genome-wide upper bound on exposure: **41/494 = 8.30% [6.18, 11.07]** of shipped gorilla families are
held together **exclusively** by edges with block ≤ 1,200 bp and coverage < 0.60 (153/2445 intra-family
edges). ⚠ That is an **upper bound, not a count** — real short-gene families (e.g. GWFAM272
LDHA/LDHB/LDHC) live in the same stratum.

### The flagship case, proven with the locus builder removed from the circuit

**GWFAM210 (gorilla): an MRPS17 hub joined to MDM2, RFTN2, GREB1, EIF3J, TRAPPC2, PIGX.**

The decisive test bypassed our pipeline entirely: rebuilt **CHM13 RefSeq curated human transcripts**
and ran the frozen mirror between them. **Five `E_r` edges form between unrelated curated transcripts:**

| pair | identity | coverage | block |
|---|---|---|---|
| GREB1 × MRPS17 | 0.7976 | 0.5415 | 1043 |
| MRPS17 × TRAPPC2 | 0.7894 | 0.5409 | 1015 |
| MRPS17 × RFTN2 | 0.7503 | 0.5392 | 1136 |
| MDM2 × MRPS17 | 0.7679 | 0.5331 | 1029 |
| MDM2[10298-11239] × MRPS17 | 0.7955 | 0.5263 | — |

Every edge uses **MRPS17[726-1695]** — ~970 bp, **54.4% of the 1,784 bp transcript**. The element is
Alu-like: MRPS17[763-1026] self-matches MRPS17[1094-1338] at id 0.765, and ~300 bp sub-units match
CEP19, PIGX, TMEM97, GREB1, RFTN2, MDM2. Probed against all 1,415 shipped gorilla reps it returns
**96 records across 57/494 = 11.54% [9.01, 14.66] of families**, while the hub's other half
(MRPS17[0-700]) and controls return exactly **1** family each.

⭐ **The gorilla node is clean by every test** — 3 exons, 4,673 bp span, aligns to human MRPS17 at
id 0.9855 covering 1.0000 of the human transcript, and the gorilla GFF puts exactly **one** gene under
it. No mis-chain, no swallowed neighbour. **The node builder cannot be blamed. The rule admitted it.**

⚠⚠ **It reproduces across species**: the same MRPS17 element drives **HSA GWFAM185**
(MRPS17 + TRAPPC2/TRAPPC2B + IVD + KCTD2/ATP5PD + LIMS1/LIMS2/LIMS3). Not a gorilla artefact.

### The external-truth false merge

**ATP1A1 × ATP4A** (human, HGNC gene_group truth: "ATPase Na+/K+ transporting" vs "ATPase H+/K+
transporting" — different groups). Curated transcript × curated transcript: **id 0.7163, coverage
0.5689, 2,117 bp aligned ⟹ edge forms.** No node involved.

⚠ The defence *"they are deep P-type ATPase paralogs, so the truth is too fine"* **does not rescue O1**:
the same panel's accepted **positive** GFPT1 × GFPT2 scores **id 0.7295, coverage 0.5353** — *lower*
coverage than the negative. **No value of τ or c orders these two correctly.**

## 3. ⚠⚠ A BUG IN A SHARED HELPER, found en route — check anything built with it

`o1_errorcensus/mkreps.py` writes exons in **ascending genomic order** and lets `bedtools getfasta -s`
reverse-complement each one individually. **Every minus-strand transcript therefore comes out with
sense-correct exons in reversed order.** ATP4A is minus-strand, so its reference rep was scrambled —
which had produced a spurious "no edge, cov 0.0443" reading and a provisional re-classification of the
case as node-construction. **That re-classification is withdrawn.** Rebuilt in true transcript order,
the node is *exonerated in the opposite direction*: NODE5 × fixed refATP4A gives **id 1.0000,
coverage 0.9745**.

⭐ **Any earlier result built with `mkreps.py` on a minus-strand gene is suspect and must be re-checked.**

## 4. The node-construction failures (47) — the definition is innocent, the catalog is not

Two recurring pathologies, both measurable at scale:

**(a) One locus cut in two.** ZNF492 window W063: identity **exactly 1.000000** over **exactly** the
1,204 bp genomic intersection of the two node intervals — the two "copies" are the same locus split.
At scale: **28/494 = 5.67% of gorilla and 15/394 = 3.81% of human families** have two copies sharing
≥ 500 bp at id ≥ 0.99; **20/494 = 4.05%** of gorilla families have two copies whose best gene is the
**same gene instance**.

**(b) Unspliced stubs and pure-repeat nodes.** ANKHD1 window W106: three single-exon unspliced
fragments of **one** gene, one of which is a 206 bp node that is **100% soft-masked repeat**. The
emitted "family" is a gene reported as a family of itself.

⚠ Honest note: (b) shares the GWFAM210 mechanism — coverage on `min(len)` lets a repeat clear `c` on a
short node. The classification split is real, but the two categories are not fully independent.

## 5. Verdict

> **O1's definition does not survive unscathed: ~30 of 105 classified bad-family cases are definitional,
> they share a single mechanism — the scale-free min-length coverage denominator admitting a ~1 kb
> dispersed repeat as ≥50% of any short node — and it reproduces in both human and gorilla with the
> locus builder removed from the circuit.**

**What survives:** identity has never been the culprit (0 failures across four substrates; the
offending edges sit at 0.75–0.80). The mechanism is **localised and named**, not diffuse. The measured
false-merge rates remain low (1.33% human/window; 6.67% human/gene-pair against external truth).

**The strongest objection still open:** no retuning fixes it. Raising `c` or `τ` cannot separate
ATP1A1 × ATP4A (negative, cov 0.5689) from GFPT1 × GFPT2 (positive, cov 0.5353) — the positive scores
*lower*. A fix has to change the **denominator or the substrate**, not the threshold. ⚠ And with the
three adversarial attacks unrun, **30 is a floor**.

---

# Appendix A — pathology (a) dissected

*Merged from `o1_error_case_census.md` on 2026-08-20. Pathology (a) is a category OF this census, so it
belongs here rather than beside it.*

## Census pathology (a) decomposes — and its largest class was misattributed

**Status 2026-08-19. Offline (T8), whole catalog, nothing through the shipped binary.**
Scratch: `/mnt/linuxdisk/home/juanfraitu/o1_gmult/patha.py`.

## 1. Why this was worth running

Today's read-based guards all died, and the reason converged: in the shipped catalog the 47
node-construction failures are dominated by **pathology (a), "one locus emitted as two"**, which is a
**coordinate/identity** signature, not a read signature. The census found it but nobody had checked
what merging the flagged pairs would do.

Four signatures, applied within each shipped family, over all 3,751 same-family copy pairs:

| signature | definition | pairs | families | census |
|---|---|---:|---:|---|
| `A_OVERLAP` | genomic intervals intersect | 35 | 35 = 7.09% | — |
| `A_ADJACENT` | same chrom, gap ≤ 10 kb | 50 | **28 = 5.67%** | **matches 28/494** |
| `A_IDENT` | rep alignment ≥ 0.99 id over ≥ 500 bp | 204 | 106 = 21.46% | broader than the census's |
| `A_SAMEGENE` | both best-overlap the **same GFF gene instance** | 29 | **20 = 4.05%** | **matches 20/494 exactly** |

Two signatures reproduce the census's numbers exactly, from an independent re-derivation.

## 2. ⭐ The largest class is NOT pathology (a)

**All 35 `A_OVERLAP` pairs are on OPPOSITE STRANDS — 35/35.** They are not one locus split in two.
They are **overlapping antisense / nested gene pairs**, a real and common arrangement, and the
`distinct_locus_reps` predicate is right to treat opposite-strand loci as distinct.

The question then becomes why they are in one *family* at all. They share the same DNA read in
opposite directions, so they align — **on the minus strand**:

| | minus-only edges | rate |
|---|---:|---:|
| the 35 overlapping antisense pairs | **33** | **0.9429** |
| all shipped edges | 56 / 2,727 | 0.0205 |

**A 46× enrichment.** These are reverse-complement artifacts, which is precisely the class the
**transcript-orientation guard** rejects. Merging them would be flatly wrong — it would fuse two
genuinely different genes. They should lose their *edge*, not their node.

### Corroboration for the orientation guard, from a new direction

The guard's evidence so far was the 74-pair FP arm (29/74 rejected, 4 lost edges of 9,032). This is
independent — it comes from **coordinates and strand**, not from the FP arm — and it shows:

* the guard's genome-wide reach is bounded: minus-only edges are **56/2,727 = 2.05%** of the catalog,
  so **that is its maximum cost**, against junction-crossing's 12.80% and the genome-anchored veto's
  3.67%. **It is the cheapest of the three guards by a factor of ~2–6.**
* **33 of those 56** are overlapping antisense pairs ⟹ **59% of what the guard removes carries no
  homology evidence.** ⚠ **Stated precisely (the earlier "provably artifactual" overstated it):** a
  gene and its antisense partner share the same DNA *by construction*, so a `-` alignment between their
  transcript-oriented reps is **entailed by the overlap** and is therefore not evidence of homology.
  That is not the same as proving the two are non-homologous — a duplicated region could in principle
  contain both. The edge is uninformative, not disproven.

## 3. What real pathology (a) costs to fix

Merging every flagged pair, judged on **family-level** outcomes only (**T7** — three prior node
changes passed on node metrics and failed end to end):

| signature | families touched | copies lost | families DISSOLVED (n < 2) |
|---|---:|---:|---:|
| `A_SAMEGENE` | 20 | 25 | **6** |
| `A_ADJACENT` | 28 | 50 | 10 |
| `A_OVERLAP` | 35 | 35 | 30 ⚠ *should not be merged at all* |
| `A_IDENT` | 106 | 148 | 75 ⚠ *too broad* |
| ANY | 133 | 201 | 86 = 17.4% of the catalog |

**`A_SAMEGENE` is the defensible detector.** It reproduces the census exactly, it is targeted (6
families dissolve, 1.21% of the catalog), and merging is the *right* operation for it: two nodes whose
best-overlapping gene is the **same gene instance** are one locus by construction.

⚠⚠ **But it reads the ANNOTATION**, and O1's standing position is *annotation as SEED, not DEFINITION*.
So `A_SAMEGENE` may be a **catalog QC flag**, never a membership condition — otherwise O1 stops being
annotation-independent, which is a larger loss than 20 families.

⚠ **`A_IDENT` must not be used.** At 106/494 = 21.46% it cannot separate a split locus from a **real
tandem duplicate at 99.9% identity** — the two present identically on identity alone. This was stated
before the numbers were computed and the numbers did not resolve it.

## 4. Net

Pathology (a) is smaller and better-behaved than "47 node-construction failures" suggests:

* **35 families** are misattributed — antisense overlap, an **edge** problem the orientation guard
  already covers, not a node problem;
* **20 families** (4.05%) are genuine same-gene-instance splits, fixable at a cost of 6 dissolved
  families, but only with the annotation;
* **28 families** (5.67%) are adjacent-locus candidates, unresolved;
* the identity-only route is unusable.

**Nothing here changes the definition.**

## 5. ⭐ CONFIRMED END-TO-END ON THE REAL BINARY (2026-08-20) — this is no longer T8

The prediction above was offline. The genome-wide catalog built with the guard **on by default** now
measures it directly:

| | families with an overlapping same-family pair | opposite-strand |
|---|---:|---:|
| 494 catalog, guard OFF | **35/494 = 0.0709** | **35/35** |
| **627 catalog, guard ON** | **4/627 = 0.0064** | 3/4 |

**A 91% reduction in the antisense-overlap class, measured by the shipped binary.** Every one of the
35 pre-guard cases was opposite-strand, which is what made the mechanism identifiable in the first
place. Together with the human negative panel (spurious E_r edges **28 → 3**,
`o1_false_positive_rules.md`) the guard's precision benefit is now measured on **two independent
substrates and the real binary**, not on the GGO FP arm it was derived from.

---

# Appendix B — the non-circular evaluation, run (2026-08-20)

Reproducer: `bench/o1_compara_noncircular.py`. HSA/CHM13 catalog (394 families, 1,220 copies) against
**Ensembl Compara paralogy**, entirely from the on-disk cache (0 live lookups).

## Why this answers the circularity objection — and what it does not

The standing objection is that seeding O1 with an annotation makes any discovery circular. That
objection is precise and it is about the **evaluation**, not the method: using annotation as a seed is
fine; scoring against truth **derived from that annotation** is not. Compara paralogy comes from
Ensembl's gene trees and is independent of the CHM13 GFF used to seed and to name nodes.

⚠ The annotation is still used here — to map a node to a gene symbol. What it is **not** used for is
deciding whether two genes are paralogues. That distinction is the whole point.

## Result

| | |
|---|---|
| copies mapped to a protein-coding gene | 987/1,220 = 0.8090 |
| same-family pairs mapping to the **same** gene (a split locus, not a paralogy claim) | 15 |
| **PRECISION** — checkable pairs Compara confirms as paralogues | **193/552 = 0.3496** |
| **RECALL** — Compara paralogue pairs the catalog puts in one family | **194/1,186 = 0.1636** |

⭐ **The headline finding is a replication.** The previously recorded 4/12 = 33% precision against
Compara now reproduces at **0.3496 on 552 pairs — 46× the scale.** It was not small-n noise.

Stratified by pair identity:

| identity | checkable | confirmed | rate |
|---|---:|---:|---:|
| ≥ 0.95 | 2 | 2 | 1.0000 |
| 0.90–0.95 | 4 | 3 | 0.7500 |
| 0.80–0.90 | 132 | 52 | 0.3939 |
| < 0.80 | 359 | 118 | 0.3287 |

The gradient runs the right way, **but ≥0.90 holds only 6 pairs in total, so nothing may be claimed
there.** 491 of 497 pairs sit below 0.90 identity.

## ⚠ Why the number is a bound on AGREEMENT, not a precision

**Compara paralogy and "multi-copy gene family" are different predicates**, and the disagreement is
concentrated exactly where they differ. Worked example, verified by hand: `ZNF684`, `ZNF26` and `ZNF84`
sit in one emitted family; each has **62, 176 and 36** Compara paralogues respectively; and **none is
listed as a paralogue of either other**. Compara's relation is tree-derived and specific; `E_r` is
sequence homology. Neither is wrong — they answer different questions.

Checked before believing the number: the cache is healthy (512/553 paralogue queries carry ≥1
paralogue, 41 empty, 0 malformed), a failed symbol lookup **skips** a pair rather than scoring it as a
miss, and the join was verified by hand on the ZNF examples.

## What to take from it

1. ✅ **A non-circular evaluation is possible and has been run.** That is the answerable form of the
   advisor's objection, and it now has a number.
2. ⚠ **Compara is the wrong external truth for this claim.** Its predicate is not O1's, so it bounds
   agreement at ~0.35 and cannot rise much regardless of how good O1 gets.
3. ⭐ **The external truths that DO match the predicate are the published expansions** — Yoo/Rhie's
   gorilla-specific MAPKBP1 / SPTBN5 / PLA2G4B, which the instrument recovers **exactly 8/9/9** at
   identity 0.973–0.983 and stable across coverage floors 0.05→0.50 — and SD catalogs. Those are the
   cards to play against the circularity objection, not Compara.
