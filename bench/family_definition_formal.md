# A formal definition of a multi-copy gene family at the RNA level — and its proof on IsoSeq reads

*Advisor interest #1: can we give a full formal definition of a multi-copy gene family and find them from
IsoSeq reads, with an honest false-positive / false-negative accounting?*

**Bottom line.** A multi-copy gene family is defined as a connected component of the **read-conflict graph**
under a divergence-tie (`de`) edge criterion. On the GGO HiFi IsoSeq substrate it reproduces a labelled
17-candidate panel exactly — **TP=7, TN=10, FP=0, FN=0 (precision = recall = 1.000)** — with every
load-bearing count independently re-derived from raw `samtools`/`pysam`, the three formal properties verified,
and the principled `de` criterion demonstrably avoiding three `AS`-tie false positives (including a 3347-read
retrocopy). The single most important caveat: the clean ledger holds for families with **≤6 cross-mapping
copies** — `minimap2`'s default secondary cap (N=5) truncates the input and fragments a real ≥12-copy NC_086018.1 (chr23)
array; this is a false negative of the *input*, not the definition, fixed by re-aligning with uncapped
secondaries.

Artifacts: `bench/family_definition_demo.py` (+ `.tsv`/`.json`), `bench/family_definition_figure.png`,
verification workflow `wf_6aa71f9e-c3d` (4 adversarial checks + synthesis). Definition predicate is the shipped
Rust criterion (`read_conflict.rs`, de-tie, commits `2e9922e..2da06c2`).

---

## The read-conflict graph and the family object

Let the substrate be a set of long reads $R$ (GGO HiFi IsoSeq, polyA-selected) aligned to a reference with
`minimap2` (primary + secondary records retained), each alignment carrying the gap-compressed divergence tag
`de:f`.

**Vertices.** $V$ = de-novo *expressed* loci — the per-transcript intervals the upstream `detect_and_assign`
pipeline emits, not annotated gene spans. (This vertex resolution is load-bearing; see P1 and Limitations.)

**Best-overlap placement.** For a read $r$ and a locus $i=(c_i,s_i,e_i)$, let $\mathrm{place}(r,i)$ be the single
alignment of $r$ on chromosome $c_i$ that maximizes overlap length with $[s_i,e_i]$ (or $\varnothing$ if none
overlaps). This is a *function*: at most one placement per (read, locus). Let $de_i(r)$ be its `de:f` divergence.

**De-tie conflict predicate.** Fix $\Delta = 0.005$, $\mathrm{DE_{max}} = 0.05$. A read $r$ *conflicts* on the
ordered pair $(i,j)$, written $\mathrm{conf}(r,i,j)$, iff both placements exist, are physically distinct
(same-locus guard: not $c_i=c_j \wedge |s_{\mathrm{place}(r,i)}-s_{\mathrm{place}(r,j)}|<200\,\text{bp}$), and
their divergences are **tied at the HiFi error floor**:
$$
|de_i(r) - de_j(r)| \le \Delta \quad\wedge\quad \max\!\big(de_i(r),\,de_j(r)\big) \le \mathrm{DE_{max}}.
$$
The read fits *both* loci comparably; `minimap2` cannot decide. This uses raw divergence `de`, **not** the
aligner's composite score `AS` — the latter folds in length and is the source of the avoided false positives.

**Edges.** With $\mathrm{MIN\_READS}=3$,
$$
(i,j)\in E \iff \big|\{\,r\in R : \mathrm{conf}(r,i,j)\,\}\big| \ge 3.
$$

**Family.** A **multi-copy gene family** is a connected component $C$ of $G=(V,E)$ with $|C|\ge 2$. This is the
*identifiability-relevant* family — the maximal set of loci among which reads are genuinely confused, exactly
the unit on which the downstream copy-assignment problem (Canzar-style conflict resolution) operates. It is
explicitly **not** a claim about evolutionary paralogy: a true paralog whose reads place uniquely (RFPL1/2/3)
is correctly excluded, and a retrocopy whose reads are decidable (EEF1A1) is correctly excluded.

## Verified formal properties

**(P1) Domain-sharer exclusion is by construction, not by threshold.** Because $\mathrm{place}(r,\cdot)$ is
single-valued, a read over a shared exon of two nested/overlapping genes has the *same physical alignment
record* selected as best overlap for both loci. Verified on GGO.bam: for CREB1~METTL21A all 192 shared reads,
and GCA~KCNH7 all 429 shared reads, the best-overlap record for both loci is coordinate-identical
($\Delta s = 0$; distinct placements $=0$). They contribute **zero** conflict pairs regardless of $\Delta$.
Domain/UTR sharing cannot create an edge — this is the failure mode of the sequence-similarity definition, here
excluded structurally.

**(P2) Exact decomposition.** Across the 7 recovered families, **0** reads have a de-tied placement pair
connecting loci in two *different* multi-locus components. The assignment problem separates exactly across
families with no shared-read information crossing component boundaries — the property the conflict graph requires.

**(P3) $\Delta=0.005$ is a data-derived error floor, not a tuned similarity threshold.** Within-family per-read
divergence has median $0.0019$, p90 $0.0211$ (the HiFi floor). Within-family conflict-read $|de_i-de_j|$
medians: MAGd2 $0.0000$, MAGd3 $0.0000$, AK6 $0.0012$, MAGd1 $0.0005$, CCDC196 $0.0037$, RABL2 $0.0061$.
Resolvable decoys sit above: CNN2-retro $0.0084$, EEF1A1-retro $0.0168$, APOBEC3 $0.0241$. A $\Delta$-sweep
reproduces ground truth for all $\Delta\in[0.003,0.006]$; first failure at $\Delta=0.007$ (CNN2), then $0.009$
(EEF1A1), then $0.017$ (APOBEC3). $\Delta=0.005$ is centered in a genuine valley with two-sided margin.
*Honesty note:* RABL2's recovery is **quorum-carried, not tie-tight** — its per-read $|de\,\text{diff}|$ median
$0.0061$ exceeds $\Delta$, and only 26% (51/195) of its conflict reads are individually $\le\Delta$; RABL2 is an
edge because 51 reads $\gg 3$ clear the bar. The "cluster near 0" picture is clean for AK6/MAGd2/MAGd3; for the
recent-paralog case it is $\mathrm{MIN\_READS}=3$ that carries it.

---

## Demonstration on GGO.bam — 17-candidate panel

The ledger reproduces **exactly**:
$$\boxed{\text{TP}=7,\quad \text{TN}=10,\quad \text{FP}=0,\quad \text{FN}=0,\quad \text{precision}=\text{recall}=1.000}$$

**Seven families found** (de-tie connected components, $|C|\ge 2$): `{RABL2A, RABL2B}`, `{AK6, LOC115934278}`,
`{CCDC196, LOC129526440}`, `{MAGd1a, MAGd1b}`, `{MAGd2a, MAGd2b}`, `{MAGd3a, MAGd3b}`, and `{MAGd0a, MAGd0b}`
(displayed as `{MAGEA9, MAGd0a, MAGd0b}` — MAGEA9 is the annotation name for the same physical locus as MAGd0b,
a redundant vertex, not a third member; see note).

### Per-family evidence (raw counts re-derived from GGO.bam)

| Candidate | Regime | Truth | de-conflict reads | AS-tie reads | Verdict |
|---|---|---:|---:|---:|---|
| RABL2 (RABL2A~RABL2B) | recent paralog, separate chrom | family | **47** / 195 | 190 | TP |
| AK6 (~LOC115934278) | high-identity cross-chrom copy | family | **24** / 24 | 13 | TP |
| CCDC196 (~LOC129526440) | high-identity cross-chrom copy | family | **24** / 27 | 27 | TP |
| MAGEA_dn0 (MAGd0a~MAGd0b) | co-located array, de-novo loci | family | **103** / 104 | 104 | TP |
| MAGEA_dn1 (MAGd1a~MAGd1b) | co-located array, de-novo loci | family | **83** / 101 | 38 | TP |
| MAGEA_dn2 (MAGd2a~MAGd2b) | co-located array, de-novo loci | family | **303** / 311 | 311 | TP |
| MAGEA_dn3 (MAGd3a~MAGd3b) | co-located array, de-novo loci | family | **75** / 75 | 75 | TP |
| APOBEC3 (D~F) | diverged paralog (resolvable) | not | 0 | 6 | TN |
| RFPL (1/2/3) | paralog, reads place uniquely | not | 0 | 0 | TN |
| EEF1A1_retro (~LOC109023808) | processed-pseudogene retrocopy | not | **0** / 3410 | **3347** | TN |
| CNN2_retro (~LOC129524764) | processed-pseudogene retrocopy | not | **1** / 127 | 121 | TN |
| ASDURF~ASNSD1 | domain-sharer (nested) | not | 0 | 0 | TN |
| CASP8~FLACC1 | domain-sharer | not | 0 | 0 | TN |
| CREB1~METTL21A | domain-sharer | not | **0** / 192 | 0 | TN |
| GCA~KCNH7 | domain-sharer | not | **0** / 429 | 0 | TN |
| GPR39~LYPD1 | domain-sharer | not | 0 | 0 | TN |
| MAGEA_annot (9/4/10/1) | annotated array genes | not | 0 | 0 | TN |

### The contrast: `de` avoids three `AS`-tie false positives

The legacy `AS`-tie criterion (composite alignment-score ratio $\ge 0.9$) forms **10** components, three spurious:

- **APOBEC3** — AS links 6 reads; de = 0 (D/F genuinely resolvable, $|de\,\text{diff}|$ median $0.024$).
- **EEF1A1 retrocopy** — AS links **3347** reads; de = 0. Reads are decidable: EEF1A1 `de≈0.001` @ MAPQ 30–48
  (primary) vs LOC109023808 `de≈0.017` @ MAPQ 0 (secondary), $|de\,\text{diff}|\approx0.016 > \Delta$. `AS`
  collapses this because scores are length-driven; `de` sees the divergence gap.
- **CNN2 retrocopy** — AS links 121 reads; de = 1 (below $\mathrm{MIN\_READS}$).

This is the headline for the advisor: the criterion is principled (raw divergence at the error floor), not a
tuned score ratio — and the principled choice is exactly what kills the 3 high-count false positives.

See `bench/family_definition_figure.png`: (A) the conflict graph with families as components and the avoided
`AS` false-positive edges drawn dashed; (B) the per-candidate ledger.

---

## What the definition does NOT claim, quantified

### FN class 1 — Secondary-cap fragmentation (GENUINE false negative; the honest example)

The definition is sound but its **input is truncated** by `minimap2`'s default placement cap (N=5 → ≤6
records/read). For any array with **>6** cross-mapping copies, each read reports only its 6 nearest placements,
so edges between farther co-members are never observed and one true family **fragments**.

Concrete, reproducible instance — gorilla chr23 (NC_086018.1; mislabeled chrY in an earlier draft — real chrY is NC_073248.2) tandem array at **22.09–22.69 Mb** (~600 kb,
~35 kb period, ≥12 copy vertices, ≥8 expressed): the capped conflict graph reports this single physical array
as **1 component of 8 loci PLUS 4 singletons {V0, V3, V10, V11}**.
- 49 reads hit the 6-vertex cap exactly (1 primary `de=0.0004` + 5 secondaries `de 0.0098–0.0200`, all inside
  the de-tie window).
- Of the 28 expressed-copy pairs, **18 (64%) are never co-observed by any single read**; the two
  highest-expressed copies V0 (129 primaries) and V1 (~267), only 14 kb apart, **share zero reads**.
- Dropped expressed members: **V0, V10, V11** (129 / 26 / 32 primaries) become singletons.
- Cap-causality confirmed: removing the 49 cap-reads further shrinks the component.

The definition is correct; **mitigation = raise `minimap2 -N/-p` (uncapped secondaries) before building $G$.**
This FN is *outside* the 17-candidate panel (no panel family exceeds 6 placements), so it does not contradict
the ledger — but it bounds the demo's generality to **≤6-copy families** unless the input is re-aligned.

> **MEASURED — the re-alignment, run (2026-06-21, `bench/array_secondary_cap_fix.py`).** Earlier this caveat was
> inferred from coverage structure, not run. It has now been run: the 737 array-region reads were re-aligned
> with `minimap2 -ax splice:hq -uf -N 50 -p 0.1`, and the de-tie conflict graph was rebuilt over the same 21
> copy vertices from the capped production BAM vs. the uncapped re-map.
> - **Cap confirmed at the read level:** 119/737 reads sit at exactly 6 records (1 primary + 5 secondaries);
>   uncapping lifts them to up to 52 alignments.
> - **Fragmentation healed:** the array's cross-mapping homology core grows from **5 → 11 copies** (components
>   touching the array 17 → 11; co-observed copy-pairs 10 → 16 / 210). The cap false negative is recovered.
> - **No over-merge (FP = 0):** a far-away expressed gene (RABL2B, 48 Mb) stays an isolated singleton in both
>   graphs; and every one of the 10 remaining uncapped singletons has **0 cross-mapping reads to the core even
>   uncapped** — i.e. it maps *uniquely* and is a genuinely resolvable copy the definition correctly keeps
>   separate, not a cap FN. So after uncapping the identifiability family is exactly the 11 cross-mapping copies.
>
> Conclusion: the secondary-cap FN is an **input artifact, fully removable by uncapped re-alignment at zero FP
> cost** — not a limit of the definition. (The earlier "12 vertices / ~35 kb period" was idealized; the real
> array is a dense homology head ~4–5 kb apart plus sparser, often uniquely-mapping, distal copies.)

### FN class 2 — Unexpressed / single-transcribed-copy families (CORRECT BY DESIGN)

A family in which only one copy is natively transcribed produces no cross-mapping → no identifiability problem →
correctly not a family (the stated scope: identifiability-relevant, not evolutionary). RFPL1/2/3 are expressed
(12/4/47 primaries) but place uniquely → 0 conflict → correct TN. Sub-caveat: a copy receiving only secondary
spillover (array V3: 0 primaries / 116 secondaries) is a confusion *destination*; isolated, it is dropped — the
same cap-fragmentation mode, not an independent error.

### FN class 3 — MIN_READS=3 strictness (NO FN on this substrate)

The only near-threshold pair is the CNN2 retrocopy at $n_{de}=1$, correctly a true negative. No real co-located
family fires below quorum; best-overlap attribution drops no real secondary on the panel.

## FP-robustness (quantified, genome-wide)

An unbiased scanner (4.78 M alignments; 230,354 reads with ≥2 placements at $de\le0.05$) gave 5128 edges / 457
components; each classical confound was adjudicated at exact gene/exon resolution:

- **Heterozygosity** — excluded by construction. Within-locus de-tied *distinct* placements = **0** across ACTB
  (16,674 reads), GAPDH, SEPTIN7, CASP8, CREB1, RABL2A, EEF1A1. Two haplotypes of one gene cannot make a
  *cross-locus* edge.
- **Domain-sharers / segmental-dup spillover** — at exact intervals every cross-family bridge collapses
  (LRPAP1-like↔DOK7-like $=0$) while within-family LRPAP1-likeA↔B $=3242$. A 336-node "mega-component" was an
  artifact of the **scanner's 50 kb bins**, not the definition — it vanishes at de-novo locus resolution.
- **Apparent unrelated-gene survivors** (SEPTIN7↔OCLN, EEF1A1↔KIAA1328) — exonic test shows both-host-exon reads
  $=0$; the reads land on a retrocopy sharing the genomic window. These are **true** retrocopy identifiability
  families mislabeled by host annotation — correct positives, not FPs.
- **rRNA/tRNA/mito** — polyA selection removes rRNA/tRNA; **0** mito↔nuclear de-tied edges ≥3.

**No genuine FP mode found.** Robustness is **conditional on tight per-locus vertices** — the shipped pipeline
uses de-novo transcript-level loci, so this is a guardrail to preserve; the 336-node binning artifact proves
coarse vertices *can* break the best-overlap surrogate.

## Residual hardcoded parameters (disclosed, not independently swept)

$\Delta=0.005$ sits in a validated valley but established only against the panel's 3 resolvable decoys (a
genome-wide decoy population could narrow it). $\mathrm{DE_{max}}=0.05$ and the 200 bp co-location guard are two
further constants; ~3% of within-family per-read $de$ exceed $\mathrm{DE_{max}}$, and neither was swept for its
own valley. Truth labels are human-assigned priors — the cross-mapping *evidence* was verified consistent with
each label, but ground-truth paralogy was not established from an external orthology source.
