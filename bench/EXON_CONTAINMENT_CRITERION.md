# Full exon containment as a family criterion — tested, does NOT beat aln_frac

**2026-07-02.** Tested the proposal: enforce **full exon containment** (every exon of one de-novo transcript homologous to an exon of the other) as a better anti-domain-share family criterion than the coverage-based `core_recip + aln_frac` rule. Script: `bench/exon_containment_criterion.py` (+ `.tsv`/`.json`/`exon_containment_mm_cache.tsv`). Held-out component-level 5-fold CV (seed 0) against `in_dna_loose` (DNA = validation only). Deterministic (minimap2 recompute == cache). Adversarial refuters: **0 flags**.

## Method

For each of 5571 core-present E_r edges: reconstruct per-exon coords from `denovo_skeletons` introns, splice the same genome-base sequence the shipped `aln_frac` uses, align the two transcripts (`minimap2 -c -x asm20`), map aligned intervals onto each transcript's exon partition, compute `exon_contain(A→B)` = fraction of A's exons ≥50% covered at ≥0.7 block-identity; **reciprocal = min(A→B, B→A)**.

## Result — the idea is structurally sound but FALSIFIED on this substrate

**1. Containment does not beat the shipped workhorse `aln_frac`.**
- Single-feature AUC: `core_recip` 0.815, **`aln_frac` 0.915**, reciprocal-exon-containment **0.900**. Containment sits *between* core and aln_frac.
- Held-out CV rule F1: coverage (`core & aln`) **0.480**; containment 0.395; combined (`core & aln & recip`) 0.442 — **strictly dominated**. `spearman(recip, aln_frac) = 0.83` → redundant.
- On the FP survivors the coverage rule keeps: recip AUC **0.567 ≈ chance**. Adds nothing.

**2. ⭐ The single-exon wall is INERT here.** Only **1 / 102,455** de-novo loci is single-exon; **0** edges are single-exon. MAGE-type intronless-CDS retrocopies still assemble as **2–3 exon transcripts (UTR introns)**, so containment is never trivially 1.0 for a genuine single-exon reason. *(Corollary: an "exclude single-exon transcripts" scope rule would remove ~nothing on this catalog.)* The operative degradation is instead **coarse-exon quantization**: few-exon loci quantize `recip` to {0, 0.5, 1}, which is exactly why recip (0.900) < the continuous `aln_frac` (0.915).

**3. ⭐ The flagship FP refutes the premise.** The **GSTM2-hub** over-merge **survives** containment (27 genuine edges, median recip **0.789**, aln 0.629) — it shares **extensive exons, not one bridging exon**. The "one big shared exon inflates coverage" mechanism does not drive any residual FP here.

**4. Cost: containment misfires on real paralogs.** Assembly fragmentation / UTR / exon-boundary noise makes real paralogs fail containment: 29/416 TP have recip = 0, 110/416 have recip < 0.5, and **39 TP have `aln_frac ≥ 0.72` but recip < 0.5**. Cutting on recip loses **~1 real TP per genuine FP removed**. Meanwhile 123/1773 genuine over-merges have recip ≥ 0.8.

## What the residual actually is (and it isn't structural)

The residual FPs are **not** single-exon-bridge or one-big-exon FPs an edge criterion can fix:
- **GSTM2 hub** — the GST domain spans *most* of the transcript across a real family → extensive shared exons.
- **MAGE** — a **family-cardinality** over-merge: 18 real paralogous loci at diploid CN 2, held by genuinely homologous edges. Edge-level criteria cannot split real homology.
- **allele-as-copy** — handled by the separate DEMOTE step, N/A to an edge criterion.

These are a **DNA-family-boundary / cardinality** problem (finer divergence than RNA resolves), not a structural-containment problem.

## Verdict

**Keep `aln_frac`.** Full exon containment is a reasonable structural idea and would earn its place on a substrate rich in *single-exon retrocopies + clean multi-exon paralogs* — but on the real gorilla RNA catalog it is lower-AUC, redundant with `aln_frac`, near-chance on the FP survivors, cannot cut the GSTM2 domain-hub, and cuts real paralogs at ~1:1. The current criterion is coverage-threshold-based but not flawed via the hypothesized mechanism; its companion `aln_frac` already is the shared-body-fraction separator, and it does the job better than exon containment.
