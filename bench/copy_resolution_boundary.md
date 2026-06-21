# The copy-resolution boundary: which multi-copy gene families can RNA resolve, and which it provably cannot

*The empirical capstone to the copy-assignment identifiability theory (`bench/copy_assignment_theory.md`).
Measured on GGO HiFi IsoSeq (FLNC, full-length mature mRNA) against the T2T reference. A clean two-tier
boundary: a resolvable majority (exonic PSVs), and a K=0 residual that is provably unresolvable per-read —
splice-site divergence is real in the reference but aligner-masked, so the entire K=0 residual co-quantifies.*

## Bottom line

Per-read copy assignment **works for the large majority** of multi-copy families and **provably cannot work**
for a small, stereotyped residual — and we can tell them apart in advance. The hard cases are not a broad
failure of RNA; they are recent inverted segmental duplications (largely X-linked cancer-testis-antigen genes)
whose entire mature transcript is identical between copies, so no spliced read carries copy information.

## The census (definitive prevalence)

Across all 1130 de-novo multi-copy families we formed 6362 co-located copy pairs and kept the **258 that are
*assignment-relevant*** — i.e. have ≥3 cross-mapping MAPQ-0 reads, so copy assignment is actually at issue
(97 families, 49,856 classified reads). Each pair is K=0 (per-read unresolvable) iff its cross-mapping reads
have `NM_A == NM_B` (the spliced read aligns equally to both copies' exons). Script + per-pair results:
`bench/copy_resolution_census.py`, `bench/copy_resolution_census.json`.

| class | pairs | fraction |
|---|---|---|
| **resolvable** (exonic PSVs, `NM_A ≠ NM_B`) | 200 | **77.5%** |
| K=0 (≥95% reads NM-identical) | 31 | 12.0% |
| K=0 strict (100% identical) | 8 | 3.1% |
| ambiguous (low support) | 19 | 7.4% |

**83.2%** of all 49,856 cross-mapping reads carry resolving evidence. K=0 is monotonically confined to the
highest-identity tail: K=0 pairs have **median 0.000% exonic divergence** vs **0.42%** for resolvable pairs.

**What the K=0 residual *is* (39 pairs / 11 families):** overwhelmingly the X-chromosome (NC_073247.2)
cancer-testis-antigen region — `MAGEA9` and the `LOC129529xxx` / `LOC1011xxxxx` tandem clusters (DNFAM1,
DNFAM14, DNFAM160, …), plus a few autosomal recent dups (`ZNF793`, scattered LOC pairs). These are recent
**inverted** segmental duplications. The big co-located families one might fear — protocadherins, keratins,
tubulins, APOBEC, APOL, BTN2A — are all in the *resolvable* majority.

## The two tiers

### Tier 1 — Resolvable majority (~77% of pairs / ~83% of reads)

The copies carry exonic PSVs; a read covering one fits its true copy better (`NM_A ≠ NM_B`). This is the regime
of **Theorem 2/3** (`copy_assignment_theory.md`): under the identifiability condition the true copies are the
unique minimum cover and are recovered in polynomial time. The shipped per-read PSV + junction assigner handles
this tier (sim K-ladder K≥2 → 100%, GGO silver-standard 100%). **This is the headline: the method works.**

### Tier 2 — K=0 residual: splice divergence is real in the reference but per-read-masked (~15-18% of pairs)

The exons are identical, so per-base PSVs give 0. The *introns* do diverge (4–60 SNPs/indels across the
duplicated block), and where that lands at a splice site it shifts a junction (e.g. MAGEA pair0's 3 bp
intronic indel makes copy A's intron 3702 nt and copy B's 3705 nt). One might hope a read using that junction
betrays its copy — but it does **not, per read**: minimap2's spliced aligner independently snaps every junction
to the nearest canonical `GT-AG` at *each* copy, so the read is fully canonical at **both** copies and, because
the indel lives inside the long intron (never in the spliced mature read, and absorbed without an exonic edit),
NM stays tied. Measured on the panel (`bench/splice_divergence_resolver.py`): per-read splice **resolution = 0
for every pair, including pair0** — the divergence is *detectable* (pair0's reads carry a different intron-length
chain at A vs B) but **non-directional** (it cannot say which copy a read came from). The earlier "pair0 ~33%"
was a reference-level count of junctions whose homologous site differs, not a per-read assignment.

So there is **no Tier-2 per-read rescue**: the K=0 residual collapses into Tier-3. (The reference-level splice
distinguishability remains relevant to copy-*model* / family-graph work — interest #3 — but not to assigning an
individual mature read.)

### Tier 3 — The K=0 core: co-quantify (no per-read assignment)

Exons identical **and** splice sites identical: the byte-identical block extends through the splice signals and
into the proximal intron, so every junction is `GT-AG` at both copies (`MAGEA pair2/pair3`-class). The 4–60
distinguishing SNPs all sit in introns/flanks a spliced FLNC read never covers; verified **0/485 junction-reads
resolvable**. Because FLNC reads are already full-length mature transcripts, **no read-length or chemistry
improvement helps** — the entire mature transcript is identical, so the information is simply not in the RNA.
This is the identifiability theorem's K=0 floor instantiated: the distinguishing-column set is empty.

For this core the correct output is **not** a forced per-read assignment but a **co-quantified ambiguity set** —
report the copies as a shared-evidence bundle under the family graph and estimate per-copy expression *ratios*
(with a DNA/copy-number prior), rather than pretending to assign individual reads. With Tier-2 shown to be
per-read-masked, the *entire* K=0 residual (not only the strict core) lands here: report a co-quantified
ambiguity set, not forced per-read assignments.

## Why this is sound (not an artifact)

The K=0 copies are genuinely distinct loci, verified against the T2T reference: distinct annotated genes
(`MAGEA9` vs `LOC129530018`), independent expression (own primary reads at each), **inverted** orientation
(B = revcomp of A — an isoform cannot be its own reverse-complement), and — decisively — each is a **bounded
44–75 kb homology block embedded in divergent unique genomic context**, the signature of a real segmental
duplication rather than an assembly false-duplication or same-locus comparison.

## The boundary as the result

```
multi-copy family with cross-mapping reads
        │
        ├─ exonic PSVs (NM_A≠NM_B)         → Tier 1  ~77%  → per-read PSV assignment (Thm 2/3)        [interest #2]
        │
        └─ K=0 (exons identical)            ~15-18%
              → splice divergence is reference-real but per-read-masked (minimap2 snaps junctions)
              → Tier 2 has NO per-read rescue → co-quantified ambiguity set (Tier 3)                  [interest #3]
```

This unifies the three advisor interests under one measured boundary: detection (#1) forms the family, PSV
assignment (#2) resolves the majority, allele-specific junctions (#3) extend resolution into part of the K=0
residual, and the same identifiability quantity (the distinguishing-column count) governs all three — with a
tiny, characterized, provably-irreducible core that should be co-quantified, not force-assigned.

## Reproducibility

- Census: `python3 bench/copy_resolution_census.py` (classifier) → `bench/copy_resolution_census.json` (258 pairs).
- Theory: `bench/copy_assignment_theory.md` (+ `bench/copy_assignment_theory_checks.py`).
- K=0 frontier attack + splice-site mechanism: workflows `wf_e36e46e4-aa1`, `wf_0c7be571-7ed`;
  `bench/resolution_improvement_bound.md`.

## Open / next

- Extend the census beyond the ≥3-read assignment-relevant subset for a fully exhaustive per-family K=0 label.
- Splice-divergence is per-read-masked by aligner junction-snapping (`bench/splice_divergence_resolver.py` — negative-result probe); the K=0 residual co-quantifies (Tier-3).
- Formalize the Tier-3 co-quantification (per-copy expression ratio under the family graph + copy-number prior).
