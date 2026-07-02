# O4 close — reference-absent copy candidates DNA-validated against the phased diploid

**2026-07-02.** The previous session left O4 (reference-absent divergent copy discovery) **mechanism-only** — it flagged 145 collapsed candidates + divergent/unmapped reps, and a *haploid* structural resolver split them ~26 copy / 38 allele / **65 undecidable**, but the haploid assembly *cannot* separate a het single-locus variant from a real extra copy, and no DNA was available. That was the biggest circularity risk in the adversarial reviews.

**What changed:** we downloaded the phased diploid **mGorGor1** (mat + pat, DNA-derived, independent of the RNA), and the RNA donor is confirmed **== mGorGor1**. So copy-vs-allele is now settleable **donor-exactly and non-circularly**: a real reference-absent **COPY** = an extra full-length locus in *both* haplotypes; an **ALLELE** = a *one-haplotype* difference at the same locus. Script: `bench/o4_diploid_validate.py` (+ `.tsv`/`.json`). Candidate consensus aligned (`minimap2 -cx asm20 -N50 -p0.1`) to `GGO.fasta` (= GCF primary, the RNA ref) + both haplotype `.mmi`; distinct full-copy loci counted at id ≥ 0.90, coverage ≥ 0.5. Deterministic (reference re-align md5 identical).

## Result — how the 145 (+6 divergent) resolve

| diploid label | count | % | meaning |
|---|---|---|---|
| **COPY** (extra full locus in **both** haps) | **0** | 0% | **no DNA-confirmed reference-absent germline copy** |
| ALLELE | 96 | 64% | one-hap difference / single-locus het |
| NOVEL | 42 | 28% | see split below |
| UNDECIDABLE | 13 | 9% | ref already ≥2 loci = mat = pat |

- **The 65 haploid-undecidable are resolved:** 47 ALLELE, 11 UNDECIDABLE, 7 NOVEL.
- **The 26 haploid-COPY calls all fall:** 16 NOVEL, 8 ALLELE, 2 UNDECIDABLE — they were SEDEF-partner *context*, not full-length extra copies. None survive.
- **ALLELE 96** = **92 plain single-locus het** (ref = mat = pat = 1; the co-segregating "PSV haplotypes" were the maternal/paternal **alleles** of a one-copy gene) + **3 genuine one-hap extra full copies** (het CNV) + 1 boundary artifact.
- **NOVEL 42** (refuter-corrected split) = **10 genuinely absent from ref + mat + pat** (novel-or-artifact) + **32 reference-present but below the 0.90 identity/coverage gate** (fragmentary/partial consensus — *not* novel, *not* reference-absent).
- **Reference-absent CONFIRMED in the diploid: 1** — `UNM_cl0`, the divergent/unmapped POC hit, absent from the primary (id 0.0) but full-length in the **maternal** haplotype at 99.7% → a real maternal-haplotype locus the primary missed (hap-specific, het-level). (CDK11B was excluded — it is present in the reference at id 0.9916, a coverage-threshold crossing, not truly reference-absent.)

## Robustness

**COPY = 0 at *every* setting** (identity 0.80/0.85/0.90, coverage 0.2–0.5, gap 100 k–500 k) — the negative is not a threshold/method-capacity artifact. **Sensitivity floor (stated honestly):** "0 COPY" = 0 copy at **≥ ~80% identity**; a genuinely divergent (<80%) both-haplotype germline copy would be scored NOVEL/genuinely-absent — the identity floor of a transcript-copy oracle. The NOVEL/ALLELE split *is* coverage-sensitive; only COPY = 0 is threshold-robust.

## Verdict — O4 is DNA-closed; the payoff is a firewall, not a discovery

**Yes as a method.** The previously-blocking ambiguity (haploid can't tell het from copy → 65 undecidable) is now **unblocked donor-exactly and non-circularly** — the independent DNA diploid, which never saw the RNA divergence signal that flagged the candidates, adjudicates every one deterministically. That closes the #1 circularity risk.

**But the scientific outcome is negative, and that is the honest, valuable headline: 0 of the 145 single-locus PSV-flagged O4 candidates (and 0 of 6 divergent) are DNA-confirmed reference-absent germline copies.** They are dominated by **heterozygosity** (64% resolve to ALLELE — the co-segregating PSVs were the two parental alleles of a one-copy gene), plus partial/artifact consensus (28% NOVEL, mostly reference-present-below-gate), and ref-already-models-it (9%). So this candidate set's single-locus collapsed-copy detector is ~100% allele/artifact at the transcript-copy level — exactly the "divergent variant could be a copy *or* an allele" caution the previous session raised, now **settled: they are alleles.**

The genuine reference-absent copies in this genome are the **family-level genomic ones `diploid_cn_oracle` already recovered** (~53 families with hap_CN > reference, non-circular), plus the **1** hap-specific locus confirmed here (`UNM_cl0`). Net: O4's copy-vs-allele ambiguity is **closed and the mechanism validated**; the result is a **precision/firewall** win (single-locus divergent candidates are alleles, not copies), not a discovery of new germline copies.
