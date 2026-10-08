# PREREG — ideal-scenario chromosome simulation

**Written 2026-09-21 before any arm is scored.** Question (user): *"simulate an ideal scenario using 1
chromosome with some multi-copy gene families but also some others that are single copy to prove that in
theory it should work ... remove unspecific phenomena like readthroughs and ensure the reads are complete."*

This establishes the **CEILING** of node construction + family definition under ideal input, and
attributes the gap to real-library phenomena. It is a positive control, not a rule search: no decision
rule is being selected, so no parameter is tuned on the outcome.

## Why this is not §6n0

§6n0 ([[project_npip_sim_ceiling]]) simulated **one family** (26 NPIP transcripts) and scored transcript
completeness: 25/26 complete with full-length jittered reads. It proved the ASSEMBLER has no defect on a
family that is entirely multi-copy. It did not test (a) a whole chromosome, (b) **single-copy genes as a
negative control**, or (c) the node/family endpoint. All three are the point here.

## Substrate

**chr16**, every annotated `gene`/`pseudogene` on it, all of their transcripts (multi-isoform, as in real
data). ⚠chr16 is heavily exposed, which is acceptable for a positive control that selects nothing.
Reference `chr16.fa`; annotation `chm13v2.0 RefSeq` (never `HSA_genomic.gff`).

**Strata, fixed from the SHIPPED chr16 homology graph before any arm runs:**
- **MULTI-COPY** = gene has >= 1 homology edge in `chr16.graph.tsv` (the shipped edge set).
- **SINGLE-COPY** = no edge. This is the negative control: these must come back as ONE node each and
  must not be drawn into a family.

## Arms (the ablation ladder)

| arm | reads | purpose |
|---|---|---|
| **A_ideal** | full-length, ends jittered +-0-30 bp, err 0.001, **no readthrough** | the ceiling |
| **A_trunc** | A_ideal + 5' truncation (`trunc_frac` 0.30) | cost of degradation |
| **A_rt** | A_ideal + readthrough molecules at the **measured 7.51%** | cost of readthrough |

The readthrough rate is not invented: 33,058 of 439,985 real primary MAPQ-60 chr16 reads (**7.51%**) hit
the exons of >= 2 distinct annotated genes at >= 25 bp. A_rt injects molecules that splice one gene's
transcript to a neighbouring same-strand gene's transcript at that rate.

⚠**Jitter is mandatory** — a sim whose reads share identical (chrom,pos,CIGAR) collapses under dedup
(§6n0). Alignment uses the shipped command: `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.
Assembly uses the shipped `--assemble-only --assembly-polish full` path, unchanged.

## Endpoint

§6u7's per-copy scorer, unchanged: universe FIXED on the arm's own annotated gene set, one-to-one
max-overlap claim, **collisions counted as misses for both genes**, a gene correct iff no collision and
>= 50% of its claimed node's exonic bp fall inside it. Reported OVERALL and **stratified by
multi-copy / single-copy**. Real-data baseline for reference: 705 / 1,329 = **0.5305**.

⚠This scores NODE and FAMILY construction, not isoform recovery — a node is judged by where it sits, not
by how many isoforms it carries.

## Bars — stated before looking

1. **A_ideal >= 0.95 overall AND >= 0.90 on the multi-copy stratum** ⟹ the method works in theory and the
   real-data gap is data, not algorithm.
2. **A_ideal < 0.80** ⟹ there IS an algorithmic defect, and the ideal input has exposed it. This is a
   real possible outcome and is not to be explained away.
3. **Single-copy stratum in A_ideal must be >= 0.95**; below that, node construction over-merges even
   without readthrough, which would contradict §6u7's attribution.
4. **A_rt is the attribution test**: if A_rt falls substantially toward the real 0.5305 while A_ideal is
   high, readthrough is confirmed as the dominant cause of the over-merge found in §6u7. If A_rt stays
   high, readthrough is NOT sufficient to explain the real gap and something else is missing.

No parameter is chosen from these numbers; the arms differ only in the stated read phenomena.
