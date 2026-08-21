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
