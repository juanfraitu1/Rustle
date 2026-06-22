# StringTie vs. copy-aware de-novo at multi-copy families (GGO) — what the copy layer adds, and what it doesn't

**Date:** 2026-06-21 · **Data:** T2T gorilla (`GGO`) IsoSeq · **Branch:** `vg/flow-capacity-apportionment`
**Reproduce:** `bench/stringtie_vs_copy_headtohead.py` (+ `copy_assign.py real`)
**Status:** revised after a 3-lens adversarial review that correctly deflated an
earlier, overstated draft. Read the "What this is NOT" section.

**Question (the advisor's):** does the copy-aware pipeline do *better* than StringTie
on GGO *because it captures copy information*? Standard isoform recall (gffcompare
FSM) can't answer it — the annotation collapses paralogs, so a correct copy split
looks redundant.

**Answer (honest, narrow):** the genuine, unique contribution is **per-copy read
attribution** — assigning each read to *which paralog copy* it came from — a
capability StringTie structurally lacks. It is **not** transcript recovery:
**in 60/63 collapse cases (95 %) StringTie already models both copies' transcripts**
(distinct intron chains), it just labels them one `gene_id` and offers no per-copy assignment.
The set of genuine tandem families where this matters is **small and must exclude
domain-sharer "families"** (the ZNF mega-family DSFAM0 accounts for 28/63 collapses
and is a false family by our own detection validation). Net: a real *capability* gain,
**bounded** by identifiability and by family-definition validity — not a headline
copy-count.

---

## Setup

- **StringTie:** `genome_st.gtf` = `stringtie -L -p 4 GGO.bam` (v3.0.1, de-novo) → 68,166
  transcripts (19,416 gene_ids). StringTie *does* read secondary alignments (verified:
  removing them changes output), but it has **no per-copy concept** — paralog copies
  whose transcripts overlap get one `gene_id` and no read-to-copy assignment.
- **Copy-aware:** the de-novo family + copy-assignment pipeline (`denovo_families*`,
  `copy_assign.py`; `bench/denovo_family_pipeline.md`).
- **Unit:** co-located multi-copy family = ≥3 *distinct* same-chrom loci within 5 Mb
  (78 such families; collapses are counted per shared StringTie gene_id).

## The over-split correction (and a real pipeline bug it exposed)

Spot-checking the top hit — **"PRNP, 5 copies"** — exposed that all five spans are
14,600,00x–14,615,51x: **one locus** split into 5 pseudo-copies (the locus union-find
merges isoforms only on an *identical* junction, so near-identical-junction variants stay
separate and then look like a multi-copy family). This was not rare: **495 of the 1,190
de-tie "families" (≈42 %) were entirely one over-split locus** (verified: ASTN, SDHAF,
PRNP, novels — all isoform-fragments of a single gene), and ~21 % of co-located "copies"
were fragments.

**Fixed at source** in `denovo_family_split.py` — an *output-level* dedup (merge family
members whose spans reciprocally overlap ≥50 %), applied **after** family detection so the
POA detection/decomposition is untouched: **1,190 → 695 family-class families** (495
over-split single-locus families dropped), MAGEA 15→**11**, **APOBEC3 recovered**,
RFPL/RABL2/KRT preserved, copies 3,636→2,362. (A first attempt merging loci *before* POA
was reverted — it perturbed detection and dropped genuine copies, MAGEA 15→6, APOBEC3
broken.) **All numbers below are on the fixed list**, where the head-to-head's own
over-split guard is a 0 % no-op — confirming the fix subsumes it.

## What the copy layer actually adds

**1. Per-copy read attribution (the real, unique capability).**
`copy_assign.py real` over 25 co-located families / 30,709 reads:
- **Unique-mapper agreement = 28,704/28,726 = 99.9 %** — this is the *silver-standard*,
  but it is on the **confident unique-mappers** (minimap2's MAPQ>0 reads, ~97 % here). It
  validates that our per-copy labels agree with minimap2 where minimap2 is already confident.
- **Hard reads:** 95 % of all reads get a confident copy assignment; these are *not*
  independently validated (no orthogonal ground truth) — reported, not claimed as proven.
- **Copy-specific junctions** resolve **+167 reads** PSVs alone cannot — e.g. **DSFAM42**
  (5 copies, 95 % MAPQ-0): **10 % → 99 % resolvable**.
StringTie produces none of this — it has no read-to-copy assignment at all.

**2. Copies StringTie merges under one gene_id (attribution, not recovery).**

| | Count | Honest reading |
|---|---|---|
| Collapse instances (≥2 distinct loci share one StringTie gene_id) | **63** | |
| └ **StringTie emits distinct intron chains at both loci** | **60/63 (95 %)** | → it *models* the copies; only the *label* is merged. **Attribution, not recovery.** |
| └ in **DSFAM0** (178-copy/18-chrom ZNF mega-family = false domain-sharer family) | **28/63** | excluded as a non-genuine family |
| Collapses in **genuine** (non-DSFAM0) families | **35** (mostly disjoint-tandem) | many still other ZNF |
| └ in genuine families **copy-assignment validated** | **12** | the defensible core |
| └ in genuine **non-ZNF** families (PDPR, GSTM, TRIM, MYH, RFPL, novels) | **~16** (≈6 validated) | cleanest, free of domain-sharer doubt |

Named, defensible wins (StringTie merges the copies' transcripts under one label; we
attribute reads per copy): **RFPL** (incl. RFPL2, a family-rescued single-read copy),
**GSTM, TRIM, MYH, PDPR, PCDHB, MAGEA/MAGEB, KRT**.

**3. Family-aware rescue (soft, overlapping, mostly thin).**
94 co-located distinct loci are **family-rescue-flagged** (POA-confirmed against a family
member — a mechanism StringTie lacks), and 55 distinct loci get **no StringTie model** at
all. **These two sets overlap heavily (53 of 55) and must NOT be summed** with the
collapses or each other. Most are single/two-read assemblies (the rescue stage is 85
single-read + 25 two-read; a handful have high support). Honest framing: a thin,
real-but-soft tail, not a headline.

## What this is NOT (corrections to the first draft)

- **Not transcript recovery.** StringTie models 60/63 collapsed copies; the gain is
  per-copy labeling/attribution.
- **Not 40–60 "recovered" copies.** That conflated overlapping win-rows and a
  domain-sharer family. The defensible collapse core is **~12 validated** (≈6 non-ZNF).
- **The 99.9 % is unique-mapper consistency** (easy reads), not proof on the hard reads
  the method exists to resolve.
- **StringTie is not primary-only.** It reads secondaries; it simply has no copy model.

## Why standard recall can't see even this

gffcompare scores against a paralog-collapsed annotation, so per-copy attribution reads
as redundant. The (separate) read-coherence recall win (+1,735 FSM) is an *assembler*
lever, orthogonal to copy-awareness.

## Q: more multimappers in GGO.bam?

GGO.bam was aligned `minimap2 -ax splice:hq -uf` (default cap `-N 5 -p 0.8`). Re-aligning
array reads uncapped (`-N 50 -p 0.1`, `arrayfix/`) turned **737 records → 10,513 (9,757
secondaries) across ~413 loci** — the multimapper graph the copy-aware pipeline uses to
*enumerate and assign* array copies (the 5→11 array-core recovery). StringTie reads
secondaries but has no per-copy model, so more of them give it no copy resolution; the
copy-aware side is the beneficiary. Caveat: 737 reads cost 18.5 GB RAM (full-genome
index) — do it **targeted at family loci**.

## Verdict for the advisor

*"On GGO, StringTie usually models the paralog copies' transcripts but merges them under
one gene_id with no read-to-copy assignment. The copy-aware pipeline's genuine, unique
contribution is **per-copy attribution** — 99.9 % agreement with confident unique-mappers,
plus copy-specific-junction resolution (DSFAM42 10→99 %). The set of genuine tandem
families where this is demonstrable is small (~12 validated collapse loci, ~6 free of
domain-sharer doubt) and bounded by identifiability **and** by which 'families' are real.
This is a capability StringTie lacks, not a transcript-recovery count."*

## Reproduce

```bash
PY=/home/juanfra/miniforge3/bin/python
$PY bench/copy_assign.py real               # -> copy_assign_real.out (resolvability)
$PY bench/stringtie_vs_copy_headtohead.py   # -> the tables above (incl. labeling-vs-recovery, DSFAM0 split)
```
