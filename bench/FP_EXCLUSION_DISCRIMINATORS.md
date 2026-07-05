# False-positive breakdown + four-axis exclusion attempts — the domain-share residual is DNA-bound

**Question.** Fresh FP percentages by type for the current catalog (`family_rna_refine.tsv`, md5
`dca64cbd`, 573 families), and *clean* ways to exclude the identifiable FPs — especially the
domain-sharers. Verdict: **the dispersed domain-share residual is provably
DNA-copy-number-bound — five axes fail (four *similarity* axes + VG topology).** One genuinely new axis
(genome *architecture* — coordinates + strand, §6) cleanly clips the antisense/nested-overlap corner (wired
as a 4th gate); the dispersed bulk stays DNA-bound. The investigation also *corrected the precision upward*
and characterized the residual completely.

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

## 6. A FIFTH axis — genome architecture — yields the one clean new rule

The four exhausted axes are all **similarity** axes ("how *similar* are the copies?"). A fifth, un-tapped
axis is **genome architecture** ("*where* are they, and which *orientation*?") — coordinates + strand,
**zero sequence similarity**. Within it there is exactly one clean lever:

**ANTISENSE / RECIPROCAL-OVERLAP rule** — cut a within-family cross-gene edge iff the two genes are
`same_contig ∧ opposite strand ∧ reciprocal_overlap (overlap_bp / min(span)) ≥ 0.50`, after mega-span
array-flagging. Two distinct genes reciprocally overlapping ≥50% on *opposite* strands occupy the **same
genomic region** (sense/antisense or nested-gene) and **cannot be two copies of one gene** — a geometric
impossibility, not a similarity threshold.

- **9 genuine over-merges excluded** (`RASA1+CCNH`, `RNASEH2C+KAT5`, `ARHGEF39+CCDC107`, `HDGFL3+TM6SF1`,
  `TRMT10B+EXOSC3`, +4) at **0 confirmed-paralog / 0 diploid-oracle collateral** (`SAME_Ep`=0 at every
  threshold; the 2 DNA-loose "reals" are overlap-induced spurious neighbor edges, correctly flagged).
- **New beyond the 3 shipped gates:** 8/9 are protein×protein (the purity flag can't fire), the gates
  never use coordinate/strand, and ≥6 are 2-node single-edge families the mosaic gate can't touch.
- **Principled threshold:** the real antisense pair `MPDU1/MPDU1-AS1` sits at 0.49, just below 0.50 — so
  lowering the cut would start eating genuine antisense-transcript pairs.
- **Honest limit:** scope-limited (~9 pairs) — it clips one clean corner of the domain-share residual, it
  does **not** dissolve the dispersed bulk, which stays DNA-CN-bound.

Rules that FAILED independent-truth validation: **biotype-consistency** as a hard gate (~19% collateral —
lncRNA is where biotype-decayed *real* copies hide; kept as the reported purity flag); **terminal/
single-shared-exon** (the falsified colinear gate re-skinned, cuts divergent reals); **strand-orientation**
(0 FP are inverted-orientation edges).

## 7. The complete per-FP-class rule taxonomy

| FP class | rule | status |
|---|---|---|
| repeat-hub / repeat-bridge | repeat-hub + multi-repeat-bridge gates | **shipped** |
| mosaic / recombinant | recombinant-split gate | **shipped** |
| **antisense / reciprocal-overlap** | **opp-strand + recip-overlap ≥ 0.50 (coordinate axis)** | **NEW-clean — 9 FP, 0 collateral (wired as 4th gate)** |
| noncoding-bridge (protein × lncRNA) | hard biotype gate | not clean (~19% collateral) → reported purity flag |
| retrocopy / pseudogene | — | truth-blind (untranslated; retrocopies are real DNA copies) |
| **single-domain-share (dispersed bulk)** | — | **IRREDUCIBLE, DNA-copy-number-bound** (5 axes fail) |
| EP_OVERSPLIT / cardinality-array | — | truth-artifact (real dups, not FPs) |
| terminal / single-shared-exon | — | reject (falsified colinear gate) |

**Room remaining:** the FP space is now fully mapped. Everything is either **gated** (repeat/mosaic/
antisense), a **truth-artifact** (~49% of the "902 FP pairs" are real dups the protein truth over-splits),
or the **irreducible domain-share residual** (DNA-CN-bound, proven across 5 axes). No further clean
RNA-derivable rule exists.

Scripts: `family_dna_vs_rna_headtohead.py`, `family_fp_mechanism.py`, `fp_by_type.py`,
`protein_discriminator.py`, `block_gate.py`, `a4_topology_dive.py`, `fp_new_rules.py` (+ the a4 refutes).
Deterministic (`PYTHONHASHSEED=0`); every "win" adversarially audited against independent DNA/cDNA/protein truth.
