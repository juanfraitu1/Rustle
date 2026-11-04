# ThesisVault methods audit, family-wise assembly — 2026-11-04

**Question (user goal):** does `~/ThesisVault` (30 papers, 14 topics) hold a method not yet applied that would move the measured
miss classes of the family-wise assembly on the r1521 panel (docs/GOAL_family_assembly_2026-11-04.md)?
**Method:** the two earlier audits (2026-09-28, 2026-10-02, both O1/O2/O3-framed) re-read; the "Method in brief" of the 13
assembly-relevant papers and the full texts of StringTie2, FLAIR, Greshnova 2024, Pal 2026, Yoo 2025 supplement grepped for
single-exon, junction, floor and support rules; each mapped onto the miss classes measured in r1532-r1536. No experiment run.

**Answer:** no vault method reaches a class that is still open. Two small rule changes are suggested by the tool papers (FLAIR's
end-density ends for single-exon models; a read-count exception to the canonical-junction gate); together they could reach at
most ~10 to 15 of the 212 misses and are not vault methods proper. Everything sequence-based (IsoCon, devider, Strainy,
Paraphase, DuploMap, Dishuck/Guitart margins) targets copies that differ in sequence, which the panel's remaining multimapper
misses do not (18 exon-identical, r1532).

## Miss classes (212 expressed members without a model) against the vault

| class (register) | members | vault methods that touch it | status |
|---|---|---|---|
| single-exon truth genes, own unspliced primaries, no model (r1534 A) | 51 | StringTie2 (single-exon transcripts under a coverage floor), FLAIR (TSS/TES from read-end density, 100 nt windows, no motif), LRGASP (ends validated by CAGE / poly(A) support) | the r1489 mono layer (PAS + no priming + 50 bp end window) reaches 15 (r1535 L1). FAM47A (381 reads, 365 ending within 50 bp, 0 primed, 2 with PAS) fails only the PAS test. **FLAIR's end-density criterion in place of the PAS motif** is the one untried rule; reach bounded by the non-primed misses (FAM47A, CENPVL1, TAS2R20, HCAR3-type, ~5-10). Histones (H4C3) have no PAS by biology and no tool in the vault models them. |
| unspliced reads inside the span (r1534 B4, 46 of 53) | 46 | same as above | same mechanism; gained 0 in L1; same candidate rule. |
| unspliced-read members (r1518 trade, r1534 B2) | 20 | StringTie2 collapses identical alignments and keeps the count (step 1) | = RUSTLE_DUP_ADMIT precedent (r1525: 68k chains at 1% CAT). Closed. |
| non-canonical junction rejected by pass 1 (r1536 K1) | 7 | StringTie2: unsupported splice sites moved to the best-supported site within 10 bp; FLAIR: sites snapped to annotation / short-read junctions within 10 nt | snapping = r1527 (FAILED, 9.4% of chains move); annotation/short reads are out of scope (genome-only rule). SPANXA2's junction sits exactly on the CAT exon boundary with TT-AG on CHM13: a **read-count exception to the canonical gate** (accept a non-canonical junction carried identically by >= k reads) is the only form; reach <= 3 members at k >= 5 (SPANXA2 53, FAM86B1 5, SIGLEC10 3). Not in any vault paper. |
| unique low-support chains, 1-6 reads (r1536 K2) | 31 | StringTie2 splice graph + network flow (junctions, not exact chains); FLAIR >= 3 reads per isoform; Pal 2026 >= 2% or >= 5 reads in >= 2 of 3 assemblies; Tomaszkiewicz replicate support | splice-graph widening = r1536 arm (FAIL: 30k chains at 3.6% CAT); replicate / pooled support = r1531 (FAIL on chains kept). Closed. |
| coordinate-duplicate chains from different ZMWs (r1536 K4) | 14 | StringTie2 step 1 (count duplicates) | = r1525 class, opt-in only. Closed. |
| exon-identical siblings (r1532) | 18 | IsoCon, devider, Strainy, Paraphase, DuploMap, Dishuck >= 1 mismatch margin, Guitart >= 10 AS margin | all need a sequence difference; the 18 have none at the transcript level. Closed (r1532 ceiling). |
| T-only / tied multimappers reachable at all | <= 10 | same, plus Canzar 2016 (rejected: no abstain), Deschamps 2020 (EM / rescue: quantification only) | graph threading r1532/r1533 FAIL; closed. |

## Papers with no bearing on assembly misses (checked)
Bailey 2002, Cechova 2019, Bosch 2007, Soto 2025 (annotation-based families, short-read expression), Rhie 2023 / Makova 2024
(ampliconic gene curation by mapping known sequences; no de novo model rule), Chaisson 2017 / pHapCompass / Clair3-RNA /
Huang 2026 (phasing, variants), Vollger 2019 SDA, Prodanov 2020/2022, Chen 2023/2025 (DNA copy phasing), Yoo 2025 supplement
(read bundles define loci as read sets with gaps <= 50 bp and compare assemblies by Jaccard: a locus definition, not an
emission rule; our tandem mega-loci of r1533 are the same construction's failure mode).

## Advisor-taste notes (Topics/Canzar_advisor.md)
StringTie2's hand-set multipliers and the 1,000-node cap are cited there as the pattern to avoid; both candidate rules above are
thresholds of that kind (k reads for a junction; an end-density fraction). If either is registered, the threshold must come from
a measured distribution (as the r1529 delta did), not be tuned on the panel.

## Verdict
Family-wise assembly on this panel is evidence-limited, not method-limited: the vault's reconstruction tools (StringTie2, FLAIR,
IsoCon) and copy-resolution tools (Paraphase, Strainy, devider, DuploMap) have each been applied or refuted in register rows
(1479, 1515, 1526, 1527, 1531, 1532-1536) or need a sequence signal the remaining misses lack. The two untried rules are listed
above with their bounded reach; neither is pre-registered here.
