# Validating the readthrough rule: which single-exon transcripts are not copies?

**Date:** 2026-07-09. **Script:** `bench/readthrough_rule_validation.py`. **Substrate:** `GGO_mm.bam`.

## The problem

A single-exon de-novo transcript spanning 30–250 kb is unspliced pre-mRNA / intronic pileup, not an mRNA.
When one is admitted as a family **copy** it corrupts the copy set (GSTM: a 30 kb single-intron transcript
spanning both GSTM5 and GSTM1 became a "copy" beside GSTM3) and makes read alignment quadratic (RFPL: a 128 kb
"transcript" hangs assignment past 400 s). Real intronless genes — retrocopies, TSPYL1, GSPT2, JUND, histones —
must survive whatever rule we adopt.

The rule proposed in `FAMILY_SPOT_CHECK.md` was: *a single-exon transcript that overlaps any spliced transcript
is the unspliced form.* **That rule does not survive validation.** Two refinements of it also fail. This
document records all four, because the failures are the argument for the one that works.

## What was tested

| rule | statement | verdict |
|---|---|---|
| **R1** | T overlaps any spliced transcript | too loose to be worth measuring — a nested intronless gene overlaps its host |
| **R2** | some **assembled** spliced transcript has an intron entirely inside T's span | **FAILS — sensitivity** |
| **R3** | any **read** carries a junction entirely inside T's span | **FAILS — specificity** |
| **R4** | T's span contains **≥ 5 distinct junctions**, each supported by ≥ 2 reads | **ADOPTED** |
| R5 | span / longest-contained-read ratio ("no molecule attests a 128 kb exon") | **FAILS — no separation** |

**R2 fails on sensitivity.** It flags the GSTM readthrough but **keeps the 128 kb RFPL giant** — the very case
that motivated the rule. The RFPL window is too sparsely expressed for a spliced transcript to be assembled, so
there is no intron to engulf. R2 depends on the assembler having already succeeded, which is exactly what fails
there. Any rule phrased over *assembled transcripts* inherits this hole; the rule must be phrased over *reads*.

**R3 fails on specificity.** Over 400 expressed annotated intronless genes, R3 drops **21.2%** of them —
TSPYL1, GSPT2, DERPC, ATXN7L3B, EPM2AIP1, TSPYL4, JUND. A handful of stray spliced reads cross essentially any
locus, so "any junction" is not evidence of anything.

**R5 fails outright.** Real intronless genes reach a span/read ratio of 7.1 (low expression → short reads),
while one artifact sits at 1.04 (a single 238 kb pre-mRNA read spans it). No separation at any cutoff.

## R4 — distinct junctions, read-level, containment-scoped

> A single-exon transcript T is the unspliced form of a locus iff **≥ 5 distinct splice junctions, each with
> ≥ 2 supporting reads, lie entirely inside T's span.**

A readthrough engulfs whole gene structures, so it contains *many distinct* junctions. An intronless gene
catches at most a stray one or two. Three properties matter:

- **Distinct**, not total: TSPYL1 has 51 junction *observations* inside its span but only **4** distinct ones.
- **Read-level**, not transcript-level: it works where the assembler failed (the RFPL hole in R2).
- **Entirely inside**: a gene nested inside another gene's intron sees the host's junctions *flanking* it, never
  within it. Nested intronless genes are therefore never penalised for their host's splicing.

### Sensitivity — every single-exon de-novo transcript in 30 regions (n = 13)

All 13 are readthroughs (spans 20–250 kb). Distinct junctions inside span:

```
 14  20.2kb    15 164.0kb    16 176.5kb    19  29.8kb (GSTM)    19 128.2kb (RFPL)
 23  24.4kb    24 160.5kb    44 156.7kb    57 136.4kb           69 248.4kb
 69 250.5kb    76 110.9kb   107 176.7kb
```

**13/13 flagged. Minimum observed: 14.**

### Specificity — expressed annotated intronless genes

Median distinct junctions = **0**. Worst observed across 120 sampled genes: **4** (GSPT2, DERPC). At the 60-gene
run reproduced by the script: worst = 3 (JUN, GPR61). **0 false positives.**

### Positive control — de-novo transcripts *at* intronless loci

The controls above are annotated gene spans, not the objects the rule actually filters. Re-running the pipeline
at four highly expressed intronless genes, it does assemble single-exon de-novo transcripts there, and the rule
keeps all of them:

| gene | de-novo transcript | exons | span | reads | distinct junctions | verdict |
|---|---|---|---|---|---|---|
| TSPYL1 | `DN_NC_073229.2_140474372_1` | 1 | 14.8 kb | 2080 | **4** | kept |
| GSPT2 | `DN_NC_073247.2_61926335_1` | 1 | 2.9 kb | 672 | **4** | kept |
| ATXN7L3B | `DN_NC_073234.2_89870620_1` | 1 | 8.8 kb | 1575 | **1** | kept |
| DERPC | (assembled spliced, 3–4 exons) | 3–4 | — | 582 | n/a | rule does not apply |

### Margin

**Controls max 4 · artifacts min 14.** Any cutoff in 5–13 separates them perfectly on this data. `MIN_DISTINCT = 5`
is chosen as "no control ever reached it" — the value is read off the data, not picked. The honest margin is
**3 junctions** (TSPYL1 at 4 is the closest control), not an order of magnitude.

## Caveats — what this does NOT establish

- **n = 13 artifacts, n = 120 + 4 controls.** Both sets are small, and all artifacts come from 30 regions of one
  testis sample. Prevalence and margin genome-wide are unmeasured.
- The `≥ 2 reads` support floor is a noise guard, not a derived quantity. At `≥ 1` read the control tail rises.
- No **retrocopy** appears among the 13 single-exon de-novo transcripts, so the rule has not been tested against
  the specific object it is most likely to destroy — a recent intronless retrocopy inside a family window. The
  positive control (TSPYL1/GSPT2/ATXN7L3B) is the closest available proxy: intronless, expressed, assembled as
  single-exon, and kept.
- DERPC is a reminder that "annotated intronless" and "assembled single-exon" are different things: the pipeline
  assembles it as a 3–4 exon transcript, so the rule never sees it.

## Recommendation

Adopt R4 as a filter on de-novo transcripts *before* family formation, so the readthrough never becomes a copy.
It is read-level, containment-scoped, needs no annotation, and has a clean empirical margin. Report the count in
the `copy_assign` log so a dropped transcript is never silent.

## Reproduce

```
python bench/readthrough_rule_validation.py --bam GGO_mm.bam \
    --artifact-gtf 'txdump/*.gtf' --controls single_exon_genes.tsv
```
`txdump/*.gtf` are `copy_assign --gtf --min-copies 99` dumps (no family survives ⇒ fast, transcripts still emitted).
`single_exon_genes.tsv` is every annotated mRNA with exactly one exon (`chrom start end strand gene`).

Related: `bench/FAMILY_SPOT_CHECK.md` (where the defect was found), `bench/artifact_audit.py`.
