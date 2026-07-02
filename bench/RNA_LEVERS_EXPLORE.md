# RNA-lever exploration — do any NEW RNA edge levers help the E_r oracle?

**Question.** The shipped RNA-only edge oracle keeps an `E_r` edge iff
`core_recip >= 0.26` AND `aln_frac >= 0.72` (`aln_frac` = shared spliced-exon-body
fraction; component-CV held-out AUC **0.915**, the dominant separator). Do any OTHER
RNA-derived levers (a) add held-out **lift** over `core_recip + aln_frac`, or (b) **crack
the irreducible residual** — real paralogs that share a long body but the diploid DNA
oracle calls DISTINCT families (MAGE sub-cluster `LOC129529978`/`986`, GSTM2 paralog
blocks)? DNA is used ONLY to LEARN the operating point and to VALIDATE/name the residual —
**never gates an edge.**

**Build.** `bench/rna_levers_explore.py` → `rna_levers_explore.tsv`, `.json`,
`rna_levers_l2_cache.tsv`. Population = the SAME **5580** core-present refined-`E_r` direct
edges (**416** TP) the shipped LEARN stage uses; K=5 whole-component folds reused verbatim
(`rna_only_edge_oracle.assign_folds`, seed=0). Deterministic: byte-identical JSON/TSV across
re-runs (`PYTHONHASHSEED=0`).

**Baseline reproduced.** `aln_frac`-alone pooled held-out AUC = **0.9152** (doc 0.915 ✓);
OOF-logistic `core+aln` = **0.9040**; `core`-alone 0.8153. (The candidate levers are scored
against the more conservative OOF-logistic `core+aln`=0.9040 baseline, so any reported lift
is if anything understated relative to the 0.915 single-feature bar.)

## Four candidate RNA levers — held-out AUC + lift over `core_recip + aln_frac`

| lever | cover | standalone AUC | base+lever AUC | lift | beats 0.915 |
|---|---|---|---|---|---|
| L1 splice / structural | 5571 | 0.603 | 0.9129 | **+0.0089** | no |
| L2 divergence / spread | 5571 | 0.905 | 0.9094 | +0.0053 | no |
| L3 expression | 5571 | 0.407 | 0.9032 | −0.0008 | no |
| L4 PSV / SUN density | 5580 | 0.575 | 0.8779 | −0.0261 | no |
| **ALL levers added** | — | — | **0.8841** | −0.0199 | no |

Strongest sub-features (single-feature pooled AUC): L2 `block_ident` **0.9104** / `total_frac`
**0.9104** (essentially the same as `aln_frac` 0.9152 — same alignment, redundant, not
additive); L1 `njunc_min` 0.6865, `intron_cos` 0.6276; L3 `log_expr_diff` inverted 0.4286;
L4 `npsv_min`/`nsun_min` inverted 0.31/0.41 at only 7.5% coverage (`l4_present`'s 0.7088 is
just a "both genes are catalogued multi-copy" membership proxy, not a density signal).

**No lever — alone or combined — beats the 0.915 baseline.** L1 is a small, real, non-decisive
precision nudge (+0.0089 to the fitted logistic ≈ within noise on 416 positives; `njunc_min`
favors multi-exon architecture over single-exon repeat bridges, consistent with the shipped
retrocopy-intronlessness filter). L3 is null (single pooled library → no co-expression axis,
only an abundance ratio). L4 is null-to-harmful (92% missing coverage and the density is
inverted; adding it drops AUC to 0.878).

## Residual-crack test (MAGE sub-cluster + GSTM2 + DNA-named roster)

98 residual edges touch the roster; **30** survive the `aln_frac >= 0.72` keep gate (12
DNA-real, 18 over-merge) = the over-merges the `aln_frac` axis alone cannot cut. **The FULL
shipped rule (`core_recip >= 0.26` AND `aln_frac >= 0.72`) already cuts 4 of those 18
over-merges with 0 real edges lost, so the true irreducible residual is 14 over-merges, not
18.** On the 30-edge keep subset:

- reference `aln_frac` AUC 0.653; reference **`core_recip` AUC 0.773** ← already in the rule,
  the separator to beat.
- best FULL-COVERAGE lever = L1 `njunc_min` (flipped AUC 0.766) / `intron_cos` 0.722; L2
  `block_ident` 0.685. **None beats `core_recip` 0.773; none cracks it.**
- the single higher-discrimination hit, L4 `npsv_min` (flipped AUC 0.939), is **LOW-COVERAGE
  (20/30), inverted-direction, and does NOT generalize** (adding L4 to the full population
  lowers AUC by 0.026) — explicitly **not a crack**.

The genuine over-merges are structurally indistinguishable from the MAGE/GSTM2 reals:
e.g. over-merge `GSTM2+LOC101129940` (core 0.79 / aln 1.02 / block_ident 0.98 / concentration
1.0) vs real `LOC101142904+LOC129530552` (core 0.38 / aln 1.00 / block_ident 1.00 /
concentration 1.0 / intron_cos 1.0 / exon_ratio 1.0) — same single-block, high-identity,
matched-architecture profile. The DNA family boundary here is **orthogonal to every RNA
homology / structure / expression / variant-density signal.** The diploid CN oracle confirms
the whole roster is genuine low-copy `multi_copy` (asm_hapCN 2–13), so these are real
sub-family splits, not artifacts.

## RNA-purity correction (retraction)

**The earlier "L2 = POA read consensus / NO genome bases / RNA-PURE" claim is RETRACTED.**
`denovo_transcripts.fa` is NOT a read/POA consensus: `denovo_assemble_gate.py` splices the
gorilla **genome reference (`GGO.fasta`)** at each locus's RNA-read-derived exon coordinates.
Verified byte-identical: transcript `DN_NC_011120.1_0_1` == `GGO.fetch(NC_011120.1, 0, 16411)`.
So L2's exon/locus **structure** is RNA-read-derived, but its nucleotide **base identities are
genome-reference bases** — identical provenance to the shipped `aln_frac`. L2 is therefore
**RNA-structure-gated, genome-base homology, NOT a base-pure lever**, and the previously
offered "purity-improving drop-in for `aln_frac`" takeaway is withdrawn (it was doubly vacuous:
`aln_frac` already uses the same genome-spliced source).

The purity guard is now **two-tier**: (1) feature-name disjointness from all DNA/protein/genome
labels (all 12 sub-features pass); (2) base-provenance — L1 (exon counts + intron lengths), L3
(read counts), L4 (read-called PSV/SUN counts) are **base-free** (structure/count only); L2 is
**genome-base**. `verify_l2_base_provenance()` asserts the byte-identity and records the genome
source rather than hiding it. There is **no base-pure sequence lever here**; a genome-base-free
sequence lever would require a true per-read POA consensus (not built).

## Honest verdict

**No new RNA lever should be added to the oracle. `aln_frac` (+ `core_recip`) is
near-complete.** No lever — alone or combined — beats the 0.915 held-out baseline
(best base+lever = L1 0.9129), and none cracks the irreducible residual (best separator is
`core_recip` 0.773, already in the rule; the only higher-disc hit, `npsv_min` flipped 0.939,
is low-coverage, inverted, and non-generalizing). L2 (divergence/spread) is a re-derivation of
`aln_frac` off the same alignment (redundant, not additive); L1 (splice/structural) is a real
but non-decisive precision nudge (+0.009); L3 (expression) is null on a single pooled library;
L4 (PSV/SUN density) is null-to-harmful. The MAGE/GSTM2 residual is **irreducible from RNA** —
it requires the DNA diploid oracle, exactly as characterized. This is a legitimate confirming
result: the shipped `core_recip + aln_frac` rule already captures the separable RNA signal.
