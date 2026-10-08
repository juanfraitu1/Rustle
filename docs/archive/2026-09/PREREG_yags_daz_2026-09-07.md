# PREREG — Y ampliconic gene families, DAZ first (2026-09-07)

**Written before the substrate was built.** md5 in `soto_mcl/yag/PREREG.md5`. Human only; never pooled
with gorilla (which annotates a single `DAZ1` on `NC_073248.2` — no family to test).

## Why this substrate
The Y ampliconic families are the case the thesis exists for. Over the human DAZ1/DAZ3 cluster
(`chrY:23,965,002-24,099,471`, A119b IsoSeq on CHM13) there are **493 primary alignments of which 358 are
MAPQ 0 and exactly one is MAPQ 60** — against a genome-wide multi-copy background where MAPQ 0 is 0.04 %
(blind-spot audit, 08-14). NPIP's hard set was MAPQ-60 reads with tied alignment scores; here the aligner does
not even reach MAPQ 60. If O2's machinery earns its keep anywhere, it is here; if it does not, that is the
result. DAZ is also palindromic (DAZ1/DAZ3 and DAZ4/DAZ2, inverted pairs) and carries internal tandem repeats,
so read-throughs between copies are expected — the §6fw guard is tested where it matters most.

## Substrate
CHM13 v2.0 **chrY** gene and pseudogene spans → all-vs-all `minimap2 -x asm20 -c -X -N 50 -p 0.1 -t 4` →
`mcl_families --min-exonic-bp 1 --merge-overlapping-loci --core-refine --sedef HSA_sedef_pairs.bed
--emit-units --emit-readthrough-units --bam winloci_data/A119b.t2t.bam --fasta chm13v2.0.fa`.
The duplication track covers chrY (1,185 pairs). Every threshold at the shipped defaults.
⚠ **Tissue caveat, stated before the run:** A119b is not testis, and the Y ampliconic families are
testis-specific. Reads over them may be mismapped from elsewhere rather than DAZ transcripts. This is a
statement about the ALIGNER's behaviour on ampliconic sequence, which is what O2 adjudicates; it is **not**
evidence about DAZ expression, and nothing here may be reported as such.

## Truth
The curated Y ampliconic symbols in CHM13 RefSeq: `DAZ1-4`, and the `TSPY*`, `RBMY*`, `CDY*`, `BPY*`, `PRY*`,
`HSFY*`, `XKRY*`, `VCY*` sets. Symbol-derived and independent of anything the method computes; its known bias
is the symbol trap (an unnamed copy is invisible), so a unit outside a truth set is reported as an unnamed
candidate, never silently as a false positive.

## Predictions
| # | prediction |
|---|---|
| **P1** | the four DAZ copies land in **one** cluster |
| **P2** | ≥ 90 % of the DAZ family's molecules are **contested** (MAPQ < 60) — this substrate is the near-tie population, not the MAPQ-0 sliver |
| **P3** | O2 assigns **≥ 20 %** of DAZ molecules to a single copy under the origin certificate, the rest abstaining. Below that, the machinery does not earn its keep here and I will say so |
| **P4** | read-through units between DAZ copies are **mostly rejected by the §6fw guard** (their flanks are duplicates by construction): **≤ 2** survive inside the DAZ cluster |
| **P5** | no molecule is assigned in violation of the tie invariant (`tie_invariant` clean in `quant.tsv`) |
| **P6** | at least one other Y ampliconic family (`TSPY`, `RBMY`, `CDY`) is recovered as a single cluster, showing DAZ is not a one-off |

## Interpretation fixed in advance
- P3 failing ⟹ report it as a limit of O2 on ampliconic sequence, not as a tuning target.
- P4 failing ⟹ the guard does not generalise from the autosomal case and must be re-examined before it ships on.
- Any claim about DAZ biology is out of scope by the tissue caveat above.

---

## AMENDMENT 1 (2026-09-07, before any result was produced)

The first build was **killed after 24 minutes of wall time and 1h44m of CPU with an empty PAF**. Cause
diagnosed, not guessed: chrY's 1,382 annotated spans include several models of **0.8–1.1 Mb**
(`LOC124908779` 41.49–42.60 Mb, `LOC124908781` 33.54–34.55 Mb, `LOC124908800` 32.01–32.84 Mb, and others),
all inside the **Yq12 heterochromatic satellite block**, and an all-vs-all of satellite against satellite is
quadratic in anchors.

**Substrate restricted to `chrY` spans starting below 28 Mb** — the euchromatic male-specific region. It keeps
every ampliconic family the pre-registration names (TSPY ≈ 6.0–6.6 Mb, RBMY ≈ 21.0–21.6, DAZ 23.97–25.72,
CDY1 ≈ 25.9–26.1) and drops the satellite models. ⚠ **No alignment parameter is changed** — `-x asm20 -c -X
-N 50 -p 0.1` and every `mcl_families` threshold stay exactly as elsewhere, so the catalogs remain comparable;
only the node set is smaller. Predictions P1–P6 stand unchanged.
