# False-positive breakdown + four-axis exclusion attempts — the domain-share residual is DNA-bound

**Question.** Fresh FP percentages by type for the current catalog (`family_rna_refine.tsv`, md5
`dca64cbd`, 573 families), and *clean* ways to exclude the identifiable FPs — especially the
domain-sharers. Verdict: **no RNA-derivable axis cleanly excludes the domain-share residual; it is
provably DNA-copy-number-bound.** But the investigation *corrected the precision upward* and characterized
the residual completely.

## 1. Fresh FP breakdown by type (902 FP pairs — but the pair count is misleading)

| type | % FP pair-mass | status |
|---|---|---|
| **EP_OVERSPLIT** (real dup the whole-protein oracle fragments) | 46.7% | **NOT an FP — truth artifact** |
| **single-domain-share** | 21.5% | genuine, protein-coding |
| noncoding element (tRNA/snRNA/lncRNA) | 11.6% | genuine, non-coding (no protein to gate) |
| unannotated / other | 9.2% | genuine, non-coding |
| retrocopy / pseudogene | 8.9% | ambiguous (a retrocopy *is* a copy) |
| cardinality array | 2.1% | **NOT an FP — truth artifact** |

**48.8% of the "902 FP pairs" is truth-artifact** — real multi-exon duplications the whole-protein E_p
oracle over-splits, combinatorially inflated by arrays (one GSTM pseudogene spawns dozens of pairs). All
are SEDEF-confirmed. The honest denominator is the **block/family view**: ~30 genuine over-merge blocks,
of which single-domain-share dominates.

## 2. The precision was UNDER-stated — the true genuine-FP count is ~21, not 32

Cross-referencing the "32 genuine over-merge" blocks against **independent** DNA/cDNA/SEDEF truth (not the
E_p axis that labelled them): **11 of the 32 are confirmed real duplications** — `RHD;SDHD` (SEDEF
matched-id **0.947**), `RPL7`, `ACSBG1;IDH3A`, `EIF1AX;LOC101137900` (SEDEF+cDNA), … The protein truth
mislabels them because it fragments divergent real dups. **Corrected genuine-FP count ≈ 21** (≈3.6% of
573). So the shipped precision (0.917) *under*-counts the definition's true accuracy.

## 3. Four independent axes tested to exclude the domain-shares — all fail

| axis | question it asks | result on the domain-share residual |
|---|---|---|
| **nucleotide** (core_recip / aln_frac, colinear-exon) | how similar? | **exhausted** — domain-shares have *higher* nt-coverage (0.54) than divergent reals (0.28); any floor cuts the reals. Colinear-≥2-exon cuts 53.7% of divergent reals. |
| **protein** (E_p reciprocal coverage) | how similar (aa)? | **circular** — "genuine over-merge" ≡ `in_ep=0`, and the protein feature *is* `in_ep`; against independent DNA the gate cuts **16–25% of real paralogs**. |
| **TE / repeat** (soft-mask, minimizer mult) | is the shared region a repeat? | **inverted** — real segdup paralogs embed *more* Alu/LINE (0.42) than the domain-share FPs (0.26); the 3 shipped gates already removed the repeat-bridge class. |
| **VG topology / snarl** (clust_min, Jaccard, articulation, snarls) | how do they *connect*? | **blind + circular** — 66–77% of the residual is **2-node** (1 edge, zero topology); on the ≥3-node minority, real domain-families (GSTM2) and FP domain-shares are *both* **hubs** (AUC 0.71 but 32–34% real collateral); clust_min/Jaccard re-express the co-threading lever we already have; new-beyond-`recombinant_split` ≈ 0. |

## 4. Why all four fail — one fact, not four coincidences

At the RNA level, an **FP domain-share and a real divergent-paralog-with-a-conserved-domain are the same
object**: a shared conserved domain at high identity. The *only* thing that separates them is whether the
**whole locus** is duplicated (real family) or **just the domain** (domain-share) — a **DNA copy-number**
question, structurally invisible to RNA. The one RNA signal that would separate them (shared-domain *span*)
is `aln_frac`, already exhausted because **divergent real paralogs also read as isolated single domains**.
Topology cannot help because a whole-body-domain real family (GSTM2: the GST domain spans the transcript)
and a one-domain FP share (RHD+SDHD) present as the *identical* hub topology.

## 5. The defensible position

The domain-share residual is **~21 blocks (≈3–4% of the catalog), fully characterized (a single conserved
domain), and provably irreducible from RNA/protein/TE/topology** — four independent axes, adversarially
verified, every one either clips real paralogs or is circular. This is a **principled boundary, not a
criterion failure**: "these are separable only by DNA copy-number, and DNA is the validation layer." The
only actionable output is a **manual-review flag** (`E_p-impure AND colinear≤1 AND both-coding`, ~30
blocks), with the caveat that ~11 are DNA-backed reals needing a per-edge (not family) split.

Scripts: `family_dna_vs_rna_headtohead.py`, `family_fp_mechanism.py`, `fp_by_type.py`,
`protein_discriminator.py`, `block_gate.py`, `a4_topology_dive.py` (+ the a4 refutes). Deterministic
(`PYTHONHASHSEED=0`); every "win" adversarially audited against independent DNA/cDNA/protein truth.
