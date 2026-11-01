# O3 — reference-absent copies: the evidence, and what it does and does not support

**2026-09-03**, updated through 2026-10-25 (registers 1408, 1416-1420). Full text: `archive/2026-11/O3_STATUS.md`. Newer O3 state: `LEDGER_O3.md` `## Current state`.

## Current state (2026-11-01)
- **Defensible sentence:** O3's RNA detectors are refuted, vacuous or unvalidated; O3 itself is UNTESTED IN THE REGIME WHERE IT WOULD OCCUR.
- **Framing (2026-10-09, Amendment 48):** O3 is DETECTION of the copy the reference lacks; placement is DNA's job. `o3_candidates` SHIPPED OPT-IN (`--candidates`); A13 deletion held-out PASS (92%, false moves 0.77%), A14 no-deletion control FAIL (35/53 false flags), fix PARKED. Real reads: recall 0.566 pat / 0.646 mat on reference-absent loci; ~4 in 5 flags are other apes' transcripts, read-share screen removes 80%. Export `bench/unmapped_rescue/o3_export.py` (KB3781: 2,559 candidates, 1,526 with read share >= 0.02).
- **How far the donor is from the reference (1408):** Jim's reference is a mosaic of his haplotypes; copy count differs in 27-35 / 494 families (5-7%), so ~30 target families on fibroblast RNA. Trio WGS: 17/516 autosomal families differ.
- **SD / exon level (1416-1420):** haplotypes differ in 44.8% of SD windows, exons of 24% of families; WGS-confirmed core 3.0% SD / 6.3% families (1419). Detectors FAIL: BPCLUST (1417, precision 0.39), per-exon certified tier (1418, 0.29). With true haplotypes 41.9% of exon gains are RNA-visible (1420): the detector, not the data, is the limit.
- **Since (O2 side, not O3 verdicts):** register 1477-1484 concern placement/O2 (`--tie-ms` opt-in, 1483); FAMASM multi-copy assembly FAIL (1479). No new O3 detector result. Open: items in `## 5. WHAT WOULD ACTUALLY SETTLE IT`.

## ⭐⭐ THE HEADLINE — AND IT IS NOT THE OBVIOUS ONE
O3 is NOT refuted; NAMED DETECTORS are. Two DNA-level positive instances exist, and no shipped binary emits a reference-absent copy at defaults. Absence of output is a consequence of the gating, not a measurement. (full: archive/2026-11/O3_STATUS.md:18)

## 1. THE THREE-WAY DECOMPOSITION — "O3" NAMES THREE SEPARABLE THINGS
| | status |
|---|---|
| specific RNA statistics | REFUTED: k+1 MEC (ARI 0.7621 / 0.2718 vs 0.80), between-copy depth (1.27x vs 109.96x variance), near-tie proxy (1.49 vs 1.56), WGS k-mer multiplicity (~1.7% shared 21-mers), `collapse_gate` (fires on EEF1A1) |
| headline zeros | VACUOUS / BLOCKED: FARDIV ∧ FARCLIP 0/915 cannot fire (FARCLIP median 0.0006 vs 0.05); unmapped route power 1/35; `absent_copy` gate 5 (0.9961 vs 0.9754, rejects 172/320); `em_copy_assign` cannot form a novel consensus; S2 mixture power 0.0588 below 0.01 divergence |
| phenomenon | DEMONSTRATED at DNA level: GOLGA6L7 (positive control, PAT-derived reference); ⭐`LOC134757045` chr14, mat 5 / pat 3, two maternal units absent from the reference, never put in front of an RNA detector |
(full: archive/2026-11/O3_STATUS.md:32)

## 2. ⛔ NO SHIPPED BINARY CAN EMIT A REFERENCE-ABSENT COPY AT DEFAULTS
`--absent-copies`, `--vg-realign`, `--collapse-gate`, `--tied-seed` all default-false; the detection was never attempted, so there is no null result. (full: archive/2026-11/O3_STATUS.md:84)

## 3. ⚠⚠ THE DO-NOT-QUOTE LIST — every one of these has been retracted or is being mis-cited
| ⛔ do not say | why |
|---|---|
| "289 candidates, 0 admitted" | only the 12-family MEC panel; 8/328 = 2.4% on 73 families (§6dc) |
| "TPR 0.7965 / FPR 0.0000 / AUC 0.8924" | RETRACTED oracle statistic (§6t) |
| "8/101 = 7.92% haplotype CNV rate" | "Do not quote 7.92%"; defensible ≤ 4/101 = 0.0396, tandem-confirmed 1/101 |
| "8/9/9 vs 5/6/8 haplotype deficit" | direction survives, magnitude does not |
| "2.78% mat-vs-pat difference rate" | UNINFORMATIVE: null 0.1512 vs 0.0278 |
| "the DNA-side 0/817 is an independent second zero" | forbidden (denominator-error pattern) |
| "the NPIP chr16:80,195,301 novel gene" | it is NPIPB15 |
| "multi-copy families do not differ between individuals" | the zero is on 99 housekeeping loci (§6dq) |
| "GOLGA6L7 is a reference-absent finding" | it is the positive control |
(full: archive/2026-11/O3_STATUS.md:94)

## 4. THE ONE WELL-POWERED ZERO, AND WHY IT DOES NOT GENERALISE
§6dq: 0 positions with ≥3 alleles across 226,615 positions, 99 loci, but they are housekeeping genes selected by the both-tissues filter (49 of 3,310 testable families, the wrong 49). (full: archive/2026-11/O3_STATUS.md:110)

## 5. WHAT WOULD ACTUALLY SETTLE IT
1. Matched-tissue IsoSeq from a second gorilla (an acquisition). 2. Put `LOC134757045` in front of the RNA detector (silent in fibroblast = UNINFORMATIVE). 3. Fix gate 5 before quoting any admission rate. (full: archive/2026-11/O3_STATUS.md:122)
