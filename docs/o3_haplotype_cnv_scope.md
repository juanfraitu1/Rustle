# Scope: how often do whole gene copies differ between one individual's two haplotypes?

**Status 2026-08-19: SCOPED, NOT RUN.** Pre-registration draft — thresholds and controls below are
declared before any result.

## The question and why it is the right one

The advisor's objection to O3 is that *"the genome is an average, so every individual differs from it
unless you compare to the individual that generated it."* For our substrate the premise is wrong on
two counts: **mGorGor1 is a haplotype-resolved assembly of ONE animal (KB3781 "Jim")**, not a mosaic
of donors, and **the fibroblast IsoSeq is that animal's own cell line**. The register already scopes
the objection as *"live only for the CROSS arm and human work."*

But the counter-proposition — *"copy number is stable; only SNPs and indels differ"* — is also false,
and our own data already refutes it. On this assembly, `_pri`/`_pat` vs `_mat`:

| probe | `_pri` / `_pat` | `_mat` |
|---|---:|---:|
| MAPKBP1 | 8 | 5 |
| PLA2G4B | 9 | 6 |
| SPTBN5 | 9 | 8 |

**One animal's two haplotypes, differing by 1–3 whole gene copies.** What is missing is the **rate**.
The existing measurement — `d_ortho` nonzero on **11/816 = 1.35% [0.75, 2.40]**, 6 up vs 5 down,
sign test p = 1.0000 — carries two stated limits: it is **tandem-only by construction** (`in_pri = 1`
for 762/817 = 93.3% of loci, so a dispersed event is structurally invisible), and the stratum spans
only **1.1224% of the genome**.

This run replaces the compartment and drops the tandem restriction.

## Design

**Probes.** 5,000 randomly sampled autosomal protein-coding genes (of 22,650 in `GGO_genomic.gff`),
**genomic spans, introns included** — the validated detection floor (P3: ≥97.3% identity over ≥5.8 kb)
was established with genomic queries. Extracted from `_pri`.

**Targets.** `mat.fa` and `pat.fa`, each restricted to its **24 chromosome-scale sequences (≥20 Mb)**.
This matters: mat has 225 sequences and pat has 24, but the 24 large ones are **98.67% of mat and
100.00% of pat**, so restricting makes completeness symmetric and discards only 48 Mb = 1.33%.

**Statistic.** `d_hap(g) = copies(g in pat) − copies(g in mat)`, counted genome-wide (not
window-restricted) at the P3 floors.

**Primary readout.** Fraction of genes with `d_hap ≠ 0`, with CI, plus the **sign test** — polymorphism
predicts symmetry, a systematic assembly deficit predicts skew.

## Controls, declared before any result

| # | control | must-pass |
|---|---|---|
| 1 | **Span-matched random intervals** — same SIZE distribution, not count-matched (count-matching alone has already killed a metric here) | the FP floor must sit well below the observed rate |
| 2 | **Known answer** — MAPKBP1 / PLA2G4B / SPTBN5 | pat 8/9/9 and mat 5/6/8, or this is not the instrument that produced the published result |
| 3 | **Self-null** — pat vs pat | `d = 0` on every probe; any nonzero is a bug |
| 4 | **Single-copy panel** — TBP, POLR2A, GTF2H1, SF1, TFRC, HMBS, PSMB6 | 1 and 1 |
| 5 | ⭐ **Probe-provenance stratification** | `_pri` is a per-chromosome MOSAIC — 16 chromosomes byte-identical to `_pat`, 9 to `_mat` — so the probe is paternal sequence on some chromosomes and maternal on others. **The rate must agree across the two strata.** If it does not, `d_hap` is measuring probe provenance, not biology |

Control 5 is the sharpest and has no analogue in the previous run, which was forced through `_pri`
as one of the two compared assemblies. Comparing mat to pat **directly** removes that entanglement
from the comparison, but not from the probe — hence the stratification.

## The main risk, stated in advance

Genome-wide counting is exactly where the previous module measured a **control artefact rate of
8.9–12.8%**, which is why it declined to treat `d_other` as its collapse statistic. If the true
polymorphism rate is ~1–2%, **the control floor could exceed the effect** and the run returns
nothing interpretable. Mitigation is the P3 floor (≥97.3% identity, ≥5.8 kb), which is what drove
`d_ortho`'s own FP floor to 0/817. **If control 1 lands above the observed rate, the correct
report is "uninformative", not a number.**

Second limit: this measures **polymorphism**, not collapse. The two are different questions and one
run cannot separate them. It answers the advisor's objection; it does not by itself advance O3's
detector.

## Cost

Two minimap2 index builds (~10 min each, ~12 GB RSS) plus ~130 Mb of query against each haplotype:
**≈2–2.5 job-hours, one at a time**, peak RSS ~14 GB of the 25 GB ceiling.

## The decision it produces

| outcome | reading |
|---|---|
| rate low (~1–2%), symmetric, control floor well below | the objection is quantitatively small; the matched design is sound and O3's negatives stand |
| rate high (>10%) | copy-level polymorphism is common ⟹ every O3 claim must be matched-individual — which we can do, and the cross/testis arm must be re-scoped |
| control floor ≥ observed rate | uninformative; report as such and do not quote a rate |
