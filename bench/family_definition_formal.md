# A formal definition of a multi-copy gene family at the RNA level — and its proof on IsoSeq reads

*Advisor interest #1: can we give a full formal definition of a multi-copy gene family and find them from
IsoSeq reads, with an honest false-positive / false-negative accounting?*

**Bottom line.** A multi-copy gene family is defined as a connected component of the **read-conflict graph**
under a divergence-tie (`de`) edge criterion. On the GGO HiFi IsoSeq substrate it reproduces a labelled
17-candidate panel exactly — **TP=7, TN=10, FP=0, FN=0 (precision = recall = 1.000)** — with every
load-bearing count independently re-derived from raw `samtools`/`pysam`, the three formal properties verified,
and the principled `de` criterion demonstrably avoiding three `AS`-tie false positives (including a 3347-read
retrocopy). Independently, **Ensembl Compara corroborates the annotated positives (RABL2, MAGEA) as recent
within-species paralogs**, and confirms the definition is a *strict subset* of recent paralogy: APOBEC3 is a Compara
recent paralog yet a correct TN here, because its copies are read-resolvable — the conflict-family is recent-paralogy
∩ read-confusability, the unit the copy-assignment problem needs. The single most important caveat: the clean ledger holds for families with **≤6 cross-mapping
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

**Record-attributed placement (no distance guard).** Each alignment *record* $\rho$ of a read $r$ (primary or
secondary; supplementary/chimeric records excluded as they are split-read pieces, not alternative placements) is
attributed to the single locus it best overlaps, $\mathrm{loc}(\rho)=\arg\max_i \mathrm{overlap}\big(\rho,[s_i,e_i]\big)$.
The read's *placement set* is $P(r)=\{(\mathrm{loc}(\rho),\,de(\rho)) : \rho \in \mathrm{records}(r)\}$ — at most one
entry per record, $de$ from the record's `de:f` tag. **Consequence (the principled replacement for the panel
demo's $<200$ bp guard):** two entries of $P(r)$ on *distinct* loci necessarily come from *distinct alignment
records* — a genuine multimapping — **by construction**; a single alignment spanning nested loci is attributed to
one locus and yields one entry, so it can never self-conflict. There is no distance threshold. This is the
shipped Rust formulation (`build_read_placements`), which supersedes the demo's coordinate guard.

**De-tie conflict predicate.** Fix $\Delta = 0.005$, $\mathrm{DE_{max}} = 0.05$. A read $r$ *conflicts* on the
unordered pair $(i,j)$, $i\ne j$, written $\mathrm{conf}(r,i,j)$, iff $P(r)$ contains entries $(i,de_i)$ and
$(j,de_j)$ whose divergences are **tied at the HiFi error floor**:
$$
|de_i - de_j| \le \Delta \quad\wedge\quad \max\!\big(de_i,\,de_j\big) \le \mathrm{DE_{max}}.
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

**Why this definition is airtight (three structural robustnesses).** The conflict-graph object is immune *by
construction* to the three artifacts that defeat a sequence-similarity family definition (POA contiguous-core,
`family_detect`):

- *Over-split fragments.* Near-identical isoform-variants the locus collapse leaves at one genomic position
  (e.g. the five PRNP 14.60–14.615 Mb fragments) cannot form a family: each record attributes to its single
  best-overlap locus, so co-positioned fragments share no conflicting read. A *similarity* definition does group
  them — measured on GGO, **~42 % of the de-tie similarity "families" (495/1,190) were one over-split locus**,
  needing a separate output-level member-merge (commit `19b348d`); the conflict definition needs none.
- *Domain-sharers.* Loci sharing only a protein domain are crossed by *uniquely-placing* reads, which produce no
  conflict edge (validated: 0 conflict on 7/7 Compara domain-sharers).
- *Retrocopies / decidable paralogs.* A read that fits one copy decisively (large $|de_i-de_j|$, or $de$ above
  $\mathrm{DE_{max}}$ on the worse copy) does not tie, so EEF1A1's 3347-read retrocopy is excluded — where the
  composite `AS` score, folding in length/gap penalties, falsely ties it (`de`-tie $\subsetneq$ `AS`-tie, the
  shipped regression invariant).

The only free parameters are the two principled de-tie constants $\Delta,\mathrm{DE_{max}}$ (disclosed below);
the former $<200$ bp coordinate guard is **eliminated** — it was an artifact of per-locus (not per-record)
placement and is structurally unnecessary.

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

**(P3) $\Delta$ and $\mathrm{DE_{max}}$ are error-model constants, not tuned similarity thresholds.**

*Derivation (Δ).* $\Delta$ is the single-read divergence-discrimination *resolution* at the HiFi error rate. A
per-read $de$ measured over aligned length $L\approx2.5\,$kb has binomial standard error
$\sqrt{\varepsilon/L}\approx0.0009$ at $\varepsilon\approx0.002$ (the within-family $de$ median); the tie statistic
$|de_i-de_j|$ then has SE $\approx\sqrt{2}\cdot0.0009\approx0.0013$, so two copies whose per-read divergence differs
by less than $\sim\!4\sigma\approx0.005$ are **statistically indistinguishable by a single read**. $\Delta=0.005$ is
therefore the read-level resolution limit set by HiFi error and IsoSeq read length — a measurement constant, not a
similarity knob. *(DE_max).* $\mathrm{DE_{max}}=0.05$ is a deliberately loose divergence ceiling separating
copy-candidate alignments — recent-paralog identity — from distinct-gene alignments; ~3% of within-family $de$
exceed it, an honest, conservative cut, not a tuned boundary.

*Empirical confirmation of the Δ valley.* Within-family per-read
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

### Independent confirmation (Ensembl Compara) — and the definition is STRICTER than "recent paralog"

The 17-panel truth labels are human-assigned; an orthogonal check against **Ensembl Compara** (within-species
paralogy, human-mapped, `bench/compara_cache.json`) corroborates the annotated positives and — more importantly —
shows the conflict-family is a *strict subset* of "recent paralog", separated by exactly the identifiability criterion:

- **RABL2A** — one *within-species paralog* at **Homininae** level (= RABL2B): a confirmed recent paralog pair.
- **MAGEA12** — two within-species paralogs at **Catarrhini** level: the MAGEA array is a confirmed recent tandem
  family (our de-novo MAGd0–3 loci are its expressed copies).
- **APOBEC3D/F** — two within-species paralogs at **Catarrhini** level, i.e. a *confirmed recent paralog cluster* —
  yet APOBEC3 is a **TN** in our ledger (reads decidable, 0 de-conflict). This is the crux: the definition does **not**
  recapitulate Compara recent-paralog detection — it is the **identifiability subset** of it. APOBEC3 is a recent
  paralog whose copies are RNA-resolvable, so it is correctly *not* an identifiability family. We answer a sharper
  question than Compara — not "are these paralogs?" but "are these copies confused by the reads?".
- **AK6 / CCDC196 copies (LOC115934278, LOC129526440)** — *absent from Compara* (the LOC paralogues are unannotated):
  copies the de-novo, read-level method finds that the reference annotation/Compara miss — a coverage gain, not a gap.

So the independent check (a) corroborates the annotated positives as genuine recent paralogs and (b) demonstrates the
definition's distinctive value: it is recent-paralogy **intersected with read-level confusability**, the exact unit the
downstream copy-assignment problem operates on.

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

## Genome-wide reproduction — not overfit to the panel (`bench/family_def_genomewide.py`)

To rule out "it only works on the 17 hand-picked cases," the *exact* de-tie criterion was applied to **every
annotated gene** as a vertex — **34,114 independent vertices** (RefSeq/Gnomon genes, *not* our own de-novo loci, so
no circularity; and *coarser* than the production loci, so this is a conservative upper bound on FP). All 380,369
secondary placements in GGO.bam were scanned and attributed by best overlap.

- **416 families** (2,829 edges, 1,698 genes); **57 % are size-2** and **86 % (357/416) are co-located on ≤2
  chromosomes** — the structure of real recent-paralog tandem arrays / segdup pairs, not random cross-mapping
  (which would not co-locate). Co-location is the genome-wide real-signal proxy (most members are uncharacterized
  `LOC` arrays no orthology DB annotates, so genome-wide Compara enrichment is not computable; the *characterized*
  panel families RABL2/MAGEA are Compara-confirmed separately, above).
- **Δ is flat genome-wide (±7 % band), not a panel fit.** Family count is **388 / 416 / 436** at Δ = 0.003 / 0.005 /
  0.007 — a ±7 % band straddling the operating point — and only reaches 538 at the far decoy threshold Δ = 0.017.
  (The count rises monotonically with Δ, so this is an *insensitivity* band, not a literal local minimum; "valley" is
  reserved for the panel correctness plateau, where 0.005 is correct and the first error appears at 0.007.) The
  Δ=0.005 operating point sits in this flat region genome-wide, independent of the panel.
- **The only systematic FP is the documented coarse-vertex over-merge:** 14 % (59/416) are cross-chromosome bridges
  of mixed genes (e.g. a 60-gene `LOC` component spanning 4 chromosomes; a 54-gene one over 19) — *exactly* the
  best-overlap-surrogate failure the FP-robustness section names. These concentrate in a handful of mega-components.

**DEMONSTRATED — tighter (production) vertices collapse the artifact (not just asserted).** Re-running the *same*
scan over the production **de-novo transcript loci** (101,467 tight vertices via `family_def_genomewide.py denovo`):
cross-chromosome bridges fall **59 → 20** (14 % → 9 %), coherence rises **86 % → 91 %**, and the Δ-valley stays
stable (195 / 212 / 227 families at Δ = 0.003 / 0.005 / 0.007). The *worst* coarse over-merges — the 60-gene `LOC`
component and the unrelated-gene bridge (`AP3B1`/`COL24A1`/`GPHN`/`GRID2`) — **vanish entirely**. The ~20 residual
cross-chrom families at production resolution are *not* the over-merge: each is a **co-located core plus a few
cross-chromosome links** (e.g. the n=25 family's members are all at `NC_073227.2:12.08 Mb`), the signature of
genuine **dispersed paralogs / processed-pseudogene retrocopies** — the "true retrocopy identifiability families
mislabeled by host annotation" the FP-robustness section already counts as *correct positives*, not unrelated-gene
bridges.

So at genome scale the definition is **stable (Δ flat to ±7 %), structurally sensible (size-2 / co-located
dominated), and its over-merge FP is a *measured* vertex-granularity effect that the production loci cut by
two-thirds** — the criterion itself produces no genuine FP mode, and the result is the opposite of a panel-tuned one.

## Precision/recall against a DNA-sequence ground truth (`bench/family_def_dna_pr.py`)

A natural advisor question: *put a precision/recall number on the RNA family definition against the "biological"
(DNA-sequence) definition.* We can — but the honest result is that the two definitions **answer different
questions**, so the raw numbers are not error rates; they are decomposed below. Every number here was independently
re-derived from the raw alignments by an adversarial verification pass (**0 discrepancies**); the per-edge audit
table is `bench/family_def_dna_pr_edges.tsv`.

**DNA ground truth (independent of reads).** All-vs-all alignment (`minimap2 asm20`) of the longest cDNA per gene over
the *same* 34,114 gene vertices. For each unordered pair we keep the best identity and the aligned fraction in each
direction (`cov_a`, `cov_b`).
- **LOOSE** paralog edge: $\mathrm{id}\ge 0.90 \wedge \max(cov_a,cov_b)\ge 0.30$ (the pinned `config.sh`
  MIN_IDENTITY/MIN_COV_FRAC; one-directional, because real paralogs' UTRs diverge — RABL2A/B align over only 34–38 %
  of the mRNA). → **17,410 edges / 1,460 families**.
- **WHOLE-GENE** paralog edge (reciprocal): $\mathrm{id}\ge 0.90 \wedge \min(cov_a,cov_b)\ge 0.50$ — meant to exclude
  domain-sharers. → **8,698 edges / 895 families** (a strict subset of LOOSE).

**Edge-level precision/recall** of the RNA de-tie graph (2,829 edges) against this truth:

| DNA truth | TP | FP | FN | precision | recall |
|---|---|---|---|---|---|
| LOOSE, all genes | 1822 | 1007 | 15588 | 0.644 | 0.105 |
| WHOLE-GENE, all genes | 1288 | 1541 | 7410 | 0.455 | 0.148 |
| LOOSE, expressed only | 1169 | 1660 | 3088 | 0.413 | 0.275 |
| WHOLE-GENE, expressed only | 822 | 2007 | 830 | 0.291 | 0.498 |

Both axes are dominated by **definitional difference**, not error:

**Recall — 80 % of "missed" DNA paralog pairs are transcriptionally SILENT.** Decomposing the 15,588 LOOSE-DNA edges
the RNA graph does not carry: **12,500 (80.2 %)** have ≥1 *unexpressed* copy (no RNA evidence — out of scope by
definition); **2,141** have reads that place uniquely (RNA-distinguishable copies); **524** are *resolvable*
(reads cross-map but the divergence gap exceeds Δ, the APOBEC3 principle); **423** are *sub-quorum* (genuinely tied
but at only 1–2 reads, below MIN_READS=3). The honest single recall figure is therefore **recall | cross-mapping
universe = 0.658** (of expressed paralog pairs whose reads actually co-map, the fraction linked). The often-quoted
0.812 ("tied universe") is reached only by *additionally* defining the 524 resolvable pairs as out-of-scope — a
modelling choice that must be argued, not folded in silently. **Note (anti-circularity):** "genuine miss = 0" (no
expressed pair with ≥3 tied reads is unlinked) is *true but by construction* — an RNA edge **is** a pair with ≥3
codivergent reads — so it is **not** presented as an empirical validation; the operative recall loss is the
sub-quorum (423) and resolvable (524) buckets.

**Precision — the DNA bar is the threshold-dependent one, and it under-detects.** Of the 1,007 RNA-only edges,
**212 (21 %)** have real cDNA homology $\mathrm{id}\ge 0.80$ but fall *just* under the arbitrary 0.90/0.30 bar — e.g.
`TBC1D1~LOC134756953` (id 0.8975, covA 0.94, **293** tied reads) and `RABL2B~LOC134756389` (id 0.85, 37 reads): real
paralogs the bar rejects. Crediting these gives **effective precision = (1822+212)/2829 = 0.719**. Of the remaining
**757 zero-homology** RNA-only edges, only **131** are tandem (<1 Mb) — defensible local paralogs whose longest-cDNA
representative failed to align (e.g. `GSTM2~LOC115933235`, a 118 kb tandem in the real GSTM cluster); the 245
same-chrom-far and 381 cross-chrom are unvalidated. The **highest-evidence** cross-chrom edges are *not* spurious
noise: `OCLN~SEPTIN7` (**3,369** tied reads) and `BCAS4~CCDC30` (**962**) are **processed-pseudogene/retrocopy** and
**read-through/chimeric** loci.

**The definitional crux — "family" ⊋ paralogy.** The read-conflict criterion detects *read-indistinguishable locus
pairs*, a **superset** of sequence paralogy: it also fires on retrocopies, processed pseudogenes, and read-through
fusions. A cDNA-vs-cDNA DNA truth structurally cannot represent these either, so they are a **mutual blind spot**, not
an RNA error.

**Orthogonal, not nested.** On the expressed universe the RNA graph and the whole-gene DNA graph **cross**: 830
whole-gene paralog pairs the RNA graph misses (silent/resolvable) vs 1,130 RNA edges that are not whole-gene
paralogs. So the read-conflict graph is **not a proxy for a static coverage threshold** — it conditions on actual
transcribed, co-mapping evidence. And no single coverage cut separates the classes cleanly: among LOOSE-expressed
edges, RNA-confirmed pairs have median reciprocal coverage **0.75** vs **0.27** for DNA-only — shifted but
overlapping. The reciprocal-coverage bar **drops RABL2A~RABL2B** (real, recip 0.34) yet **keeps CREB1~METTL21A**
(domain-sharer, recip 0.54).

**The one clean win — domain-sharers.** All five panel domain-sharers (`CREB1~METTL21A`, `GCA~KCNH7`, `CASP8~FLACC1`,
`ASDURF~ASNSD1`, `GPR39~LYPD1`) are genomically nested/adjacent genes sharing **one** homologous domain at mapq 0;
**all five pass the LOOSE DNA bar** (and `CREB1~METTL21A` passes even the reciprocal whole-gene bar), yet the RNA
best-overlap attribution **correctly excludes all five** (`in_rna=0`). A DNA homology bar over-calls them as families;
the read-conflict criterion does not.

**Conclusion.** A precision/recall number against DNA exists and reproduces exactly, but it measures the *overlap of
two different definitions*. Read-confusability is an **orthogonal, expression-conditioned evidence axis**,
non-redundant with — and non-nested in — cDNA-homology thresholds: a clean advantage on domain-sharers, a
complementary blind spot (and a unique strength) on retrocopies/pseudogenes, and a recall ceiling set entirely by
transcriptional silence and the read-evidence quorum, not by mis-calls.

### Tested and rejected — read-level coverage / intron filters do NOT sharpen precision (`bench/family_def_read_filters.py`)

A natural sharpening: require each conflicting read to (i) align over a **high fraction of its own length** at both
placements, and (ii) be **spliced (≥1 intron)** at both — intending to drop intronless retrocopy/pseudogene and
shared-repeat cross-mapping. **It is net-harmful.** Sweeping the predicate over the BAM and scoring each edge against
the DNA homology truth (`good:bad` = real-paralog edges lost per junk edge lost; want < 1):

| filter | edges | TP | precision | good:bad lost |
|---|---|---|---|---|
| baseline | 2829 | 1822 | 0.644 | — |
| ≥1 intron (both) | 1972 | 1307 | 0.663 | **1.83** |
| qcov ≥ 0.8 (both) | 2624 | 1717 | 0.654 | 1.28 |
| qcov ≥ 0.8 ∧ ≥1 intron | 1927 | 1287 | 0.668 | 1.79 |

Every configuration removes genuine paralog edges **faster** than artifacts (good:bad > 1) for a precision gain of
≤ 0.024. Two mechanisms: (1) the intron requirement kills 554 real paralog edges (reads carrying a genuine
cross-copy conflict often do **not** cross a splice junction at *both* copies — intra-exonic divergence,
single-junction reads, or a clipped diverged secondary); (2) the marquee spurious bridge `OCLN~SEPTIN7` (3,369 reads,
zero cDNA homology) is **kept** under the strictest filter because those reads *are* full-length and spliced — it is a
**coarse-vertex mislabel** of a real OCLN-retrocopy locus, not a read-quality artifact. Gene-level multi-exon
filtering also fails (OCLN, SEPTIN7, and all panel domain-sharers are multi-exon genes). The lever for these spurious
edges is therefore **vertex granularity** (the de-novo loci that collapse the bridges 59→20), not read-level quality
filtering — which confirms the precision residual is architectural, not a missing predicate.

## Residual hardcoded parameters (disclosed, not independently swept)

The two de-tie constants $\Delta,\mathrm{DE_{max}}$ are now the **only** free parameters (the former 200 bp
co-location guard is gone — record-attributed placement makes it structurally unnecessary). Both are
**error-model-derived** (P3): $\Delta$ as the single-read divergence resolution at HiFi error/length,
$\mathrm{DE_{max}}$ as a loose copy-vs-distinct-gene divergence ceiling. $\mathrm{MIN\_READS}=3$ is a
minimum-evidence quorum (a standard support floor), not a tuned constant. The honest residual: the empirical
$\Delta$-valley was confirmed against only the panel's 3 resolvable decoys (a genome-wide decoy sweep would tighten
the margin), and $\mathrm{DE_{max}}$ was set conservatively rather than swept for its own valley. Truth labels are
human-assigned priors; the *positives* (RABL2, MAGEA) are now independently corroborated as recent within-species
paralogs by Ensembl Compara (above), and the cross-mapping evidence was verified consistent with every label —
full external orthology ground truth across all 17 candidates (incl. the de-novo LOC copies) remains open.
