# The ThesisVault's methods, assessed for O2 (2026-10-20)

**Goal (user, 2026-10-20):** "look into the @ThesisVault/ for pertinent methods that could be used for O2 or add them to the
negative results if not applicable".

**O2, as the user defines it** (`feedback_o2_thesis_core_2026-10-20`): for a long read whose primary and secondary flags are
a coin toss between near-identical copies, find its true copy from PSVs, BAM tags, statistics or other means, assign or
abstain (identical copies: the set is the answer), and use the result to improve assemblies that ignore secondaries.

**Source:** all 24 paper notes in `ThesisVault/Papers/` (their "Method in brief" sections and O2 paragraphs), checked
against this repo's negative-results register. The two methods with no register row (rescue / EM; facility location) were
tested exploratorily today with `bench/o3_graph/o2prior.py` (`personal/window/prior/eval.json`); nothing here is a
pre-registered verdict.

## Verdicts

| paper (vault note) | method, as it would apply to O2 | status for O2 | evidence |
|---|---|---|---|
| Vollger 2019 SDA; Chaisson 2017 | PSVs, read × PSV matrix, nearest copy by Hamming distance; correlation clustering of reads | **Pertinent, in use.** Reference PSV columns: tested, fail on fresh data. PSVs observed in the sample's own confidently placed reads: the corrector, the only vault-lineage method with held-out gains | r1348 (Pareto set), r1350 (substitution count) fail; corrector r1343 H1 130 fixes / 3 breaks, 40 families better / 1 worse; r1346 408 / 40; r1347 82 / 0 |
| Deschamps 2020 (review): set / clustering; uncertainty reporting | report the ambiguous group instead of a copy | Tested: the set holds the truth, but a PSV-free null does as well; the decided / coin-toss certificate equals an AS null | r1348, r1350 |
| Deschamps 2020: rescue, EM, MMR (uniquely-mapped-read ratio) | give a tied read to the copy with the most uniquely placed reads | **Not applicable (exploratory today)** | see "Statistics" below; row 1351 |
| Canzar 2016: facility location | reads may only go to "open" copies (parsimony) | **Not applicable (exploratory today)** | see below; row 1352 |
| Dishuck 2025 (NPIP): best hit beats next by >= 1 mismatch on 169 haplotypes; Guitart 2024 (TBC1D3): primary beats all by >= 10 AS, else drop | a margin rule on mismatches or AS | Tested in our forms: margin on the reference = substitution / NM count (fails); population haplotypes = population PSVs (nothing portable); dropping the ambiguous = abstention (cleans one copy, loses the other's reads) | r1350, r1273, r1349; row 1353 |
| Huang 2026 longcallR; Hosseini 2025 pHapCompass; Shaw 2025 devider | phasing / haplotype reconstruction from read co-occurrence | Not applicable: phasing is outside O2 (user, 10-20); their per-read scoring against known copies = the Pareto / likelihood set already tested; devider refuted as an O3 arm | r1348; devider arm REFUTED on all six clauses (branch `devider-arm`, 6994b228, `docs/DEVIDER_ARM_*_2026-10-0x.md`); row 1354 |
| Zheng 2025 Clair3-RNA; longcallR's editing class | separate variants, RNA editing and artefacts | Not applicable: the allele-blind rules fail on genomic donor alleles, not editing (18 / 26 gorilla breaks confirmed on his haplotype; S2 break bases shared by >= 30% of reads) | r1350 audit; row 1355 |
| Sahlin 2018 IsoCon; Tomaszkiewicz 2023 | reference-free clustering + consensus, reads to the closest candidate | Not applicable to O2 (no reference copies; cannot separate identical copies); used for O3 | r1195-r1229; row 1356 |
| Pal 2026 | per-read variant fingerprint against one representative copy | Not applicable as is (never resolved to a copy); its idea, read the sample's variants, is the corrector | row 1357 |
| Kovaka 2019 StringTie2; Tang 2020 FLAIR; Pardo-Palacios 2024 LRGASP | assemblers and their benchmark | Not applicable: baselines. FLAIR's MAPQ >= 1 realignment to its isoforms re-introduces the coin toss | r1344, r1349; row 1358 |
| Greshnova 2024 (own); Makova 2024; Rhie 2023; Yoo 2025; Bailey 2002; Bosch 2007; Cechova 2019; Soto 2025 | no read-to-copy assignment (catalogs, assemblies, SD detection, family definition) | Not applicable to O2 (O1 / O3 papers) | row 1359 |

## Statistics for the reads the bases cannot decide (exploratory, `o2prior.py`)

Reads O2's set leaves with more than one copy (window 30), truth determined. **Rescue**: the kept copy with the most
uniquely placed reads (primary, MAPQ >= 20, not tied; a fragment's own source left out). **Parsimony**: assign only when
exactly one kept copy has any uniquely placed read.

| substrate | reads | minimap2 right | rescue right / decided (minimap2 on those) | parsimony right / decided (minimap2 on those) |
|---|---|---|---|---|
| S1 one genome, round 27 (unique reads) | 1,838 | 997 | 619 / 989 (539) | 111 / 229 (104) |
| HG002 Kinnex on CHM13 | 6,526 | 4,962 | 2,865 / 6,508 (4,959) | 264 / 375 (144) |
| A119b fragments, seed 2 | 1,265 | 966 | 808 / 1,251 (963) | 52 / 54 (37) |
| gorilla held-out | 1,038 | 1,019 | 1,027 / 1,034 (1,015) | 129 / 130 (129) |
| *equal-AS ties only:* S1 | 1,669 | 831 | 539 / 861 (413) | 103 / 217 (92) |
| HG002 | 184 | 68 | 103 / 182 (67) | 45 / 76 (20) |
| A119b fragments | 396 | 184 | 254 / 385 (184) | 26 / 27 (11) |

- On reads the set keeps ambiguous, rescue is **harmful** on real data (HG002: 2,094 fewer right): minimap2's primary is
  right on 76% of them, using evidence the substitution set ignores.
- On true equal-AS ties, where minimap2 is a coin toss, rescue is right 57-66% against minimap2's 37-48%: some signal,
  far from provenance. Its precision does not rise cleanly with the support ratio (S1 0.54 / 0.71 / 0.73 / 0.57 at ratio
  < 2 / >= 2 / >= 4 / >= 8; HG002 no trend; A119b fragments up to 0.90 on 63 reads).
- Parsimony fails on copies that are expressed but never uniquely placed (S1: 103 / 217 right).
- These methods are quantification devices (expected counts), not per-read provenance.

## What is left for O2

- The only vault-lineage method with held-out gains on real reads is the **sample-aware corrector** (SDA / Chaisson's
  idea, PSVs read from the sample's own confidently placed reads, no phasing). Whether it is inside O2's scope is the
  user's decision (asked 2026-10-20).
- For reads identical in the sample too (DAZ), the set is the answer; no vault method recovers provenance there.
- Not in the vault and not tried: none of the papers models the read's *structure* against copy-specific gene models
  (ends, junction sets) as evidence for exact ties; that would need a copy-specific annotation difference to exist.
