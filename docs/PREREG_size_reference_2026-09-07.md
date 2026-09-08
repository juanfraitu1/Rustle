# PREREG — the SIZE reference for bipartite matching, and the exonic-core bias test (2026-09-07)

**Written before both runs.** md5 in `soto_mcl/PREREG_size.md5`.

## Part A — stop scoring size against Soto's intervals
**User: "do not prefer Soto's size since it might be truncated, instead use the real size from the
annotation."** Soto's members are 98 %-identity SD intervals (§6fz, register 741); an SD interval can stop
where the duplication stops, not where the gene does, so a predicted unit that spans the whole gene is scored
as OVER-EXTENDED against a truncated truth. A gene-preferred bed already exists and **differs from the SD
intervals on 222 of 362 members**, with length ratios from 0.10× to 141×.

**Change:** the size reference becomes the **annotated gene span** — for each Soto member, the CHM13 RefSeq
gene or pseudogene with the greatest overlap, taken at its full annotated extent. Members with no annotated
gene (11 of 362) are **excluded from the size metric and reported separately**, never silently dropped.
The membership metrics (sensitivity, specificity) are unchanged and still use Soto's own intervals.

| # | prediction |
|---|---|
| **A1** | the in-band 0.5–2× fraction is **higher** against gene spans than against Soto's SD intervals |
| **A2** | the median size ratio moves **towards 1.0** |
| **A3** | the over-extended tail (ratio ≥ 2) **shrinks** — that is the class the SD truncation manufactures |

## Part B — does an exonic core remove the pseudogene preference?
Register 743: the core rule drops pseudogenes at **0.026** against **0.074** for protein-coding genes, i.e. it
RETAINS pseudogenes at ~3× the rate, because a pseudogene's span is mostly duplicated sequence and so clears
`core ≥ span/2` trivially. Candidate corrective: score the core in **exonic** bases — a member is kept when its
shared segment covers half its EXON length (or half the family's median exonic core).

| # | prediction |
|---|---|
| **B1** | the drop-rate ratio (pseudogene / protein-coding) moves **towards 1.0** from today's 0.35× |
| **B2** | protein-coding drop rate does **not** rise above 0.10 — the corrective must not cost genes |
| **B3** | `transcribed_pseudogene` (110 of 168 pseudogene members) remains **less** affected than non-transcribed pseudogene, since it is transcribed and therefore has exons in the core |

B1 failing ⟹ the exonic core is not the lever and the duplicon-vs-gene gap needs a different one; that will be
reported, not tuned.
