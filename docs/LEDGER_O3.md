# Ledger: O3 (copies the reference lacks)

Merged 2026-10-24 from 51 dated docs, oldest first; each kept verbatim under its own heading (its headings one
level down). The originals are in `docs/archive/2026-10/merged/` (git history has them at their old paths).

## Contents

- O3_CANDIDATES_AMENDMENT15_2026-10-05.md
- O3_MATERNAL_REFERENCE_2026-10-08.md
- PREREG_o3_maternal_reference_2026-10-08.md
- PREREG_unmapped_rescue_2026-10-08.md
- UNMAPPED_RESCUE_2026-10-08.md
- PHASED_TOY_2026-10-12.md
- PREREG_phased_toy_2026-10-12.md
- PREREG_missing_copy_posctl_2026-10-14.md
- PREREG_o3_candidates_foreign_check_2026-10-14.md
- PREREG_o3_fewread_flags_2026-10-14.md
- PREREG_o3_structure_projection_2026-10-14.md
- PREREG_o3_family_best_2026-10-15.md
- PREREG_o3_graph_real_2026-10-15.md
- PREREG_o3_graph_round4_2026-10-15.md
- PREREG_o3_graph_sim_2026-10-15.md
- PREREG_o3_graph_union2_2026-10-15.md
- PREREG_o3_graph_union_2026-10-15.md
- PREREG_o3_psv_clusters_v2_2026-10-15.md
- PREREG_o3_realdata_2026-10-15.md
- PREREG_o3_unexplained_2026-10-15.md
- PREREG_o3_graph_round10_2026-10-16.md
- PREREG_o3_graph_round11_2026-10-16.md
- PREREG_o3_graph_round12_2026-10-16.md
- PREREG_o3_graph_round5_2026-10-16.md
- PREREG_o3_graph_round6_2026-10-16.md
- PREREG_o3_graph_round7_2026-10-16.md
- PREREG_o3_graph_round8_2026-10-16.md
- PREREG_o3_graph_round9_2026-10-16.md
- PREREG_o3_graph_round13_2026-10-17.md
- PREREG_o3_graph_round14_2026-10-17.md
- PREREG_o3_graph_round15_2026-10-17.md
- PREREG_o3_graph_round16_2026-10-17.md
- PREREG_o3_graph_round17_2026-10-17.md
- PREREG_o3_graph_round18_2026-10-18.md
- PREREG_o3_graph_round19_2026-10-18.md
- PREREG_o3_graph_round20_2026-10-18.md
- PREREG_o3_graph_round21_2026-10-19.md
- PREREG_o3_graph_round22_2026-10-19.md
- PREREG_o3_npip_candidates_realign_2026-10-21.md
- PREREG_o3_npip_maternal_copy_2026-10-21.md
- O3_FINGERPRINT_SURVEY_2026-10-22.md
- O3_SCREEN_SPIKE_2026-10-22.md
- LANDING_LOCI_NEARBEST_2026-10-23.md
- O3_FIRST_PILOT_2026-10-23.md
- O3_ONE_COPY_SIM_2026-10-23.md
- O3_PORT_RESULTS_2026-10-23.md
- PREREG_o3_port_2026-10-23.md
- PREREG_o3_port_v6_heldout_2026-10-23.md
- PREREG_o3_putative_clip_p10_2026-10-23.md
- PREREG_o3_putative_heldout_2026-10-23.md
- PREREG_o3_two_tier_p4_2026-10-23.md


---

## O3_CANDIDATES_AMENDMENT15_2026-10-05.md

*(original: `docs/archive/2026-10/merged/O3_CANDIDATES_AMENDMENT15_2026-10-05.md`)*

## O3 candidates, Amendment 15/15b — the consensus-defect fix and its re-runs, 2026-10-05

Prereg: `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendments 15/15b (written 2026-10-03, before any re-run).
Defect diagnosis: `docs/archive/2026-10/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md` (byte-identical reproduction, 33 clusters of 5 families).
Code: commit `a13b817f` (the correction + regression tests); harness commit `209a5fd7` (`ACC=a15h`, env-overridable work dirs).
This doc records the implementation and every registered re-run; nothing here re-tunes a registered rule.

### The correction (as registered, no new constants)

1. **(a) cs normalisation** (`normalize_cs_ops`, before the vote): an insertion that begins with the template bases its
   adjacent skip removes (`+X·E ~|X|`, minimap2 2.30 `splice:hq`'s placement of an exon the template lacks — 107 of 134
   insertion-plus-skip pairs in the defective clusters) is shortened by them and the skip becomes matches; chained skips
   repeat.
2. **(b) one copy per identical long insertion per consensus**: carriers aggregated across columns, placed at the column
   with the most carriers (ties to the earliest); counts stay column-local.
3. **(c) the majority test**: a ≥ 20 bp insertion enters only with ≥ 3 carriers **and** `2 x count >=` the members covering
   the column — the rule the < 20 bp insertions, substitutions and deletions already obey.

`cargo test --release`: 1,112 passed, 0 failed, including the three required regression tests (both cs patterns with full
and chained skips; one identical long insertion at two columns) and the two unit tests reworked from the pre-15 rule.
Offline prediction check on the sound A13 controls: GWFAM100:c1 4,630 → 4,502 bp (−128, the predicted minority exon),
exactly the loss range the amendment measured offline (128–216 bp).

### Re-runs and verdicts (same bars, nothing re-tuned)

All runs with the corrected binary (`a13b817f`, sha1 6fc286a0), the registered substrates, no `RUSTLE_CACHE_DIR` (every
minimap2 call computed). Work dirs `a13a15/`, `a14a15/`, `a15h/` — the registered `a13/`, `a14/` products untouched.

#### A13 again (the 53-family deletion held-out) — A13-1/2/3 all PASS

- Stage: 53 families, 363 clusters, 77 candidates, **49 flagged** (vs 82 flagged in the uncorrected run), 1,802.9 s = **30.0 min** (bar ≤ 40).
- **A13-1 PASSES**: D right **7,882** ≥ 0.80 × C = 5,842.4 (107.9 % of C; identical to the uncorrected run's 7,882 — the
  correction lost no correct placements); false moves 322/41,727 = **0.77 %** ≤ 5 %.
- **A13-2 PASSES**: the unions keep **99.33 %** (pooled) / 99.43 % (isolated) of the components' reads (bar ≥ 95 %).
- The flags' composition changes exactly as the diagnosis said it should: D-derived **30/30 kept**, survivor-derived
  **46 → 13** (the 13 kept carry the divergent-but-genuine flags: their best-hit identity×coverage was already high in
  the uncorrected run), "elsewhere" 6 → 6. Deleted copies with ≥ 1 D-derived flag: 25/53 (unchanged); the cause table is
  unchanged in shape (23 no-D-reads-in-net / 2 no-cluster / 3 linked).

#### A14 again (the no-deletion control) — C2' HOLDS, **C1' FAILS**

- Stage: 53 families, 415 clusters, 40 candidates, **17 flagged in 10 families** (uncorrected: 54 class-c flags).
- Classes: 1 a_haplotype_only (GWFAM175), 2 b_allele, 14 c_unmatched.
- **C1' FAILS**: families with ≥ 1 class b/c/pri flag = **10/53 = 18.9 %** (bar ≤ 8 = 15.1 %). The 16 false flags are the
  defect doc's *other* category, not the fixed mechanism: median 11.5 reads, mismatch-dominated best hits
  (median whole-length d 0.036, median id×cov 0.974) — few-read or paralog-mixed consensus, conservative by construction.
- **C2' HOLDS**: false moves 317/59,013 = **0.54 %** ≤ 5 % (295 out of 317 from GWFAM175 alone).

#### Held-out H (Amendment 10's read set, nothing deleted; the 30 disjoint families) — H3 HOLDS, **H1 FAILS**

New `ACC=a15h` harness mode (commit 209a5fd7): `refabsent/R0.bam` (32,219 reads), W.copies restricted to the 30 disjoint
families (66 copies; GWFAM175 is dev-overlapping and correctly absent), classification by Amendment 9's rule.
- Stage: 30 families, 208 clusters, 18 candidates, **11 flagged in 6 families**; 507.7 s = 8.5 min.
- Classes: 1 a_haplotype_only, 2 b_allele, 8 c_unmatched (identical at δ/2, δ, 2δ).
- **H1 FAILS**: 6/30 families with ≥ 1 class b/c/pri flag (bar ≤ 4/30).
- **H3 HOLDS**: false moves 73/27,757 = **0.26 %** ≤ 5 %.

#### Decision rule (Amendment 15, registered): the default-on flip requires A13-1/2/3 AND C1'/C2' AND H1/H3

C1' and H1 fail → **the stage STAYS OPT-IN** (`--candidates`, ruling R14 stands); the corrected stage replaces the
uncorrected one as the opt-in arm — the duplication defect is gone (below), and the remaining false flags are a named,
different phenomenon (few-read consensus) for a future amendment, not this one.

### Reported beside (as registered)

- **Identity×coverage of every flagged union** (best genome hit, by source): A13 corrected — D median 0.9993, the 13
  survivor-derived 0.91–0.99, elsewhere 0.9988. The uncorrected run's low-id×cov survivor flags (0.445–0.89, incl.
  GWFAM37:2 0.445) are exactly the ones the correction removed; everything kept was ≥ 0.9128.
- **Duplication signature** (test (b) of the defect doc: off-diagonal self-alignment ≥ 200 bp, `minimap2 -x asm20 -X -c`):
  A13 22/82 → 14/49; A14 10/56 → 5/17. Restricted to unions below 0.99 identity×coverage (the defect's regime): A14
  6 → 3, and the 3 remaining have documented non-defect causes (GWFAM335_0 is the 204-kb genome-alignment split;
  GWFAM175_0/GWFAM331_0 are the few-read flags counted in C1'). Genuine tandem repeats survive by design — the gene
  families are segmental duplications.
- **Whole-BAM cost of the corrected stage (R23 repeated)**: one 50-family batch of the 378-family table on the full
  fibroblast BAM (23 GB), batch 0 of `wholebam/batches.txt`: **STOPPED at the 570 s limit, elapsed 631.15 s, peak RSS
  10.67 GB** (`a14a15/wholebam/wreport.out`) — byte-for-byte the uncorrected run's numbers (631.18 s / 10.69 GB,
  `a14/wholebam/wreport.out`): same pass-A/B reads, the majority test's cost is invisible next to BAM I/O + minimap2.
  Ruling R23 stands: a whole-BAM run needs smaller batches (or the second machine); the stage is not whole-BAM-ready
  at either binary.
- Stage wall times: A13 30.0 min (uncorrected 25.9), A14 27.6 min, H 8.5 min — the majority test costs ~15 %, within the
  registered 2× IsoCon budget (A13-3 bar).

### Files

| where | what |
|---|---|
| `a13a15/` | the corrected A13 run (all steps, incl. arm M, comparator, keep, report) |
| `a14a15/` | the corrected A14 run (classify, arm C, creport) + `wholebam/` R23 probe |
| `a15h/` | the held-out H run (classify, arm C) |
| `/tmp/dupsig.*.names` | the per-union off-diagonal sets behind the signature counts (session scratch) |

### Corrections (2026-10-06)

Found by an adversarial re-check of the Locus Anatomy page against the run products; each recomputed here from the products. The verdicts above (A13-1/2/3 PASS, C1' and H1 FAIL) do not change.

1. **A13-1, D right.** The uncorrected run's value is **7,898** (108.1 % of C) with 323/41,727 false moves (`a13/comparator.out`, `a13/score.out`), not 7,882. The correction cost 16 correct placements (-0.2 %) and one false move. "Identical to the uncorrected run's 7,882 — the correction lost no correct placements" (A13 section, first bullet) is wrong, and so is "unchanged by the fix" in register row 1247.
2. **A13-3, time.** The registered measure is the sum of the five batches' `/usr/bin/time` Elapsed: corrected **2,002.1 s = 33.4 min** (205.88 + 330.65 + 513.44 + 557.72 + 394.39), uncorrected **1,555.8 s = 25.9 min**: the majority test costs **+29 %**, not "~15 %". The 30.0 min above is the stage's own monotonic clock (1,802.9 s), which runs about 11 % behind Elapsed. The bar (40 min) is met on either clock.
3. **GWFAM175_0 and GWFAM331_0 are not few-read flags.** GWFAM175_0 has 183 reads (it is the class-a true copy, GWFAM175_B0) in both the A13 and A14 corrected runs; GWFAM331_0 has 143 reads in A13 and 129 in A14 (`cand.candidates.tsv`). The "few-read flags counted in C1'" sentence of the duplication-signature bullet holds for neither.


---

## O3_MATERNAL_REFERENCE_2026-10-08.md

*(original: `docs/archive/2026-10/merged/O3_MATERNAL_REFERENCE_2026-10-08.md`)*

## What happens to the reads of a copy the reference haplotype lacks, and what a reference-free step recovers (2026-10-08)

Study of `docs/PREREG_o3_maternal_reference_2026-10-08.md` (body, Amendment 1, Amendment 2); plan `docs/superpowers/plans/2026-10-08-o3-maternal-reference.md`; code `bench/o3_maternal/` (86 unit tests);
work dir `/mnt/linuxdisk/tmp/o3_mat/`; artifact https://claude.ai/artifact/R1QoJa9GLfxnt51haKG6k8 (private). Substrate: KB3781 (male) fibroblast Iso-Seq against its own maternal (`mat`) and paternal (`pat`) assemblies.
Every number is copied from `W/<ref>/fate/fate.tsv`, `side.json`, `W/<ref>/score/*.json`, `W/compare_*.json` and `W/unm.txt`. A fresh reviewer recomputed them from the BAMs and PAFs; the interpretation below was rewritten after that review.
Per prereg Q1 and S7, nothing is pooled across loci.

### 1. Results

**What the "missing" copies are.** 19 of the 20 large judged copies (7 with the mother's genome as reference, 13 with the father's) have a relative on the reference genome within the allele cutoff: nearest-relative identity 0.9911 to 1.0000 (table below). They are absent by the registered synteny rule and mostly present by sequence. GWFAM175_B0 is the only large copy whose nearest relative on the reference is farther than the cutoff (identity 0.9126), so it is the only test of whether a missing copy hides.

**Q1, where the reads of a missing copy go (per copy; mother as reference, then father).**
- No judged copy has an unmapped fate (0 reads unmapped; 1 and 8 reads partial in the two directions), but this is observed, not guaranteed: every read of these sets had already mapped to the combined primary assembly (section 3).
- GWFAM175_B0 (304 reads): absorbed (285 on another locus, 17 on the nearest relative, 2 tied) at median divergence 0.0662 against 0.0007 on the father's genome, MAPQ 60 on both: a divergence pile, not a hole.
- The other 6 judged copies on the mother's genome: REFUTED at 4 (GWFAM208_B0, 214_B0, 64_B17, 227_B0; the first three have identical or near-identical counterparts, identity 1.0000, 0.9995, 0.9999, de about 0.001 on both genomes), SUPPORTED at GWFAM205_B0 (81 of 82 reads tied; the reads' primary sits on a mat copy at 0.9913 with a second equidistant one at 0.9935), between bars at GWFAM64_B11 (3 of 21 tied).
- On the father's genome: REFUTED at 11 of 13 judged copies, PARTLY GWFAM382_B0, SUPPORTED GWFAM390_B1 (20 of 20 tied, two equidistant relatives at 0.9919 and 0.9942). All 13 have relatives at 0.9927 or closer.
- The registered verdict function uses fractions only (its clause "with de ≈ divergence" was not implemented), so REFUTED here includes copies whose reads sit on an identical copy at divergence near zero.
- Smaller expressed copies beyond the cutoff (reported, not judged): on `mat`, GWFAM175_B1 (7 reads: 1 tied, 6 absorbed elsewhere, de 0.0719) and GWFAM175_B2 (4 reads, absorbed); on `pat`, GWFAM26_B2 (17 reads, identity 0.9672: 9 tied, 8 absorbed, de 0.0102), GWFAM26_B3 (15 reads, 0.9898: 5 tied), GWFAM26_B6 (8 reads, 0.9862: 7 tied) and GWFAM26_B1 (4 reads, 0.9879: none tied). The tied-heavy ones are the closest to the advisor's mechanism; their n is 8 to 17 reads.
- Control, all other reads: 6.4% tied on `mat` (1,910 of 29,807), 5.0% on `pat` (1,461 of 29,148). A further 3,996 (`mat`) and 4,560 (`pat`) of the 35,094 reads, 11% and 13%, are tied or unplaced on the truth haplotype and are excluded from both the loci and the control.

**Unmapped reads exist, in one place.** Of the 959 reads unmapped on the combined primary assembly, 132 map on the mother's genome and 3 on the father's. 125 of the 132 lie in one 71-kb stretch of mat chromosome 12 (95.96 to 96.03 Mb, median MAPQ 60, median read length 3,607 bp; none of the 378 family copies is there): an expressed sequence the mother's genome has and the primary assembly and the father's genome do not, whose reads are unmapped without it. This is the unmapped fate, in the one read set that was not pre-filtered by the primary alignment, and it is a single locus.

**LRPAP1.** By the registered rule the only LRPAP1 locus absent from the mother's genome is the chrY copy (sex control, 1 read). p12 is not absent: it lifts as class `T?` to a mat locus 97.2% identical (Amendment 2). The chr12 cluster has 4 loci on `pat` and 3 on `mat`, and the registered lifts are many-to-one there. All 83 p12 reads are TIED with p14 even on `pat` (4 bp in 2,446); on `mat` 66 are absorbed on the p14 ortholog (identity 0.9959), 16 tied, median divergence 0.0025 against 0.0008 on `pat`. Separately, the 9 flagged IsoCon transcripts of the LRPAP1 family match pat chr12 22.55 Mb (c01, the highest-expressed chr12 copy) at identity 0.998 to 1.000, while their best mat match is 0.987 to 0.990 (mat chr3 and chr16 copies): c01's transcripts have no mat sequence within the allele cutoff, although the registered lift calls c01 present (it sends it to mat 34.29 Mb with 222 mismatches in 20 kb). So the extra chr12 locus may be c01, not p12; this study did not test it.

**Q2, recovery from the mother's genome alone.** IsoCon recovers GWFAM175_B0 (R1 PASS) and flags nothing within the cutoff (R2 PASS, 13 copies). The registered R3 FAILS (6 of 24 families without a missing copy, 25%, bar 20%), but 4 of those 6 families (GWFAM309, GWFAM318, GWFAM70, LRPAP1) are sequences found on the paternal genome at identity ≥ 0.998 whose best maternal match is 0.980 to 0.990 and that no truth locus lists, because the truth table (`bonly`) only lists loci on chromosomes the primary assembly took from the reference haplotype. Only GWFAM169 (a partial consensus at 0.95 on both genomes) and GWFAM382 (0.92 on both) match nothing. Both readings are given; the registered verdict is FAIL. The in-house stage (with the Amendment 15/15b fix) flags 1 candidate and recovers nothing (R1 FAIL, R2 PASS, R3 PASS 1 of 24). With the father's genome as reference (held-out direction) R1 is not testable (no judged copy beyond the cutoff), R2 passes for both arms (26 copies), R3 passes (IsoCon 2 of 21, in-house 0 of 21), and no candidate falls in a truth locus; flagged IsoCon candidates in GWFAM162, GWFAM256 and GWFAM70 are again sequences present on the mother's genome at ≥ 0.998 with the best paternal match below 0.99.

**Side by side.** On the genome that has the copy the reads sit on it at median divergence within the allele cutoff (≤ 0.00958) for 20 of the 21 large copies; the exception is GWFAM382_B0, whose reads are 5 to 6% divergent from both genomes (section 6, item 7). The chain's view of GWFAM175_B0: in the mother's-genome run IsoCon has 1 transcript of the copy, kept as a new copy and recovered by a candidate with 4 transcripts; in the father's-genome run the same transcript matches the reference and is not flagged.

### 2. Fate tables (LARGE = at least 20 reads; per locus)

#### Reference = mother's genome (`mat`); truth = `pat`

| locus | kind | reads | unmapped + partial | tied | on nearest relative | elsewhere | nearest-relative identity | 1 − identity | de (reference) | de (other haplotype) | verdict (fractions only) | registered prediction (REFUTED) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GWFAM208_B0 | catalog | 551 | 0 | 1 | 550 | 0 | 1.0000 | 0.0000 | 0.0008 | 0.0008 | REFUTED | held |
| GWFAM175_B0 | catalog | 304 | 0 | 2 | 17 | 285 | 0.9126 | 0.0874 | 0.0662 | 0.0007 | REFUTED | held |
| GWFAM214_B0 | catalog | 122 | 0 | 0 | 122 | 0 | 0.9995 | 0.0005 | 0.0010 | 0.0010 | REFUTED | held |
| LRPAP1_p12 | lrpap1_desc | 83 | 1 | 16 | 66 | 0 | 0.9959 | 0.0041 | 0.0025 | 0.0008 | descriptive |  |
| GWFAM205_B0 | catalog | 82 | 1 | 81 | 0 | 0 | 0.9935 | 0.0065 | – | 0.0006 | SUPPORTED | FAILED |
| GWFAM64_B17 | catalog | 31 | 0 | 0 | 31 | 0 | 0.9999 | 0.0001 | 0.0015 | 0.0015 | REFUTED | held |
| GWFAM227_B0 | catalog | 21 | 0 | 0 | 21 | 0 | 0.9911 | 0.0089 | 0.0057 | 0.0013 | REFUTED | held |
| GWFAM64_B11 | catalog | 21 | 0 | 3 | 18 | 0 | 0.9980 | 0.0020 | 0.0008 | 0.0024 | between bars | FAILED |

Control, all other reads (29,807): 193 partial, 1,910 tied, 27,704 mapped untied, median de 0.0016.

#### Reference = father's genome (`pat`); truth = `mat`

| locus | kind | reads | unmapped + partial | tied | on nearest relative | elsewhere | nearest-relative identity | 1 − identity | de (reference) | de (other haplotype) | verdict (fractions only) | registered prediction (REFUTED) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GWFAM309_B1 | catalog | 482 | 2 | 0 | 480 | 0 | 0.9991 | 0.0009 | 0.0010 | 0.0011 | REFUTED | held |
| GWFAM70_B2 | catalog | 195 | 1 | 0 | 194 | 0 | 0.9999 | 0.0001 | 0.0007 | 0.0007 | REFUTED | held |
| GWFAM402_B0 | catalog | 100 | 2 | 1 | 97 | 0 | 0.9991 | 0.0009 | 0.0018 | 0.0019 | REFUTED | held |
| GWFAM70_B4 | catalog | 97 | 0 | 0 | 97 | 0 | 0.9927 | 0.0073 | 0.0013 | 0.0090 | REFUTED | held |
| GWFAM70_B3 | catalog | 65 | 0 | 0 | 65 | 0 | 0.9998 | 0.0002 | 0.0012 | 0.0012 | REFUTED | held |
| GWFAM390_B0 | catalog | 59 | 0 | 1 | 58 | 0 | 0.9996 | 0.0004 | 0.0010 | 0.0010 | REFUTED | held |
| GWFAM169_B1 | catalog | 46 | 0 | 0 | 46 | 0 | 1.0000 | 0.0000 | 0.0010 | 0.0010 | REFUTED | held |
| GWFAM162_B0 | catalog | 41 | 0 | 0 | 41 | 0 | 0.9993 | 0.0007 | 0.0008 | 0.0011 | REFUTED | held |
| GWFAM309_B4 | catalog | 35 | 0 | 0 | 35 | 0 | 1.0000 | 0.0000 | 0.0009 | 0.0009 | REFUTED | held |
| GWFAM70_B1 | catalog | 29 | 0 | 0 | 29 | 0 | 0.9967 | 0.0033 | 0.0043 | 0.0026 | REFUTED | held |
| GWFAM169_B0 | catalog | 26 | 0 | 0 | 1 | 25 | 0.9992 | 0.0008 | 0.0008 | 0.0012 | REFUTED | held |
| GWFAM382_B0 | catalog | 22 | 3 | 2 | 17 | 0 | 0.9948 | 0.0052 | 0.0565 | 0.0570 | PARTLY | FAILED |
| GWFAM390_B1 | catalog | 20 | 0 | 20 | 0 | 0 | 0.9942 | 0.0058 | – | 0.0013 | SUPPORTED | FAILED |

Control, all other reads (29,148): 189 partial, 1,461 tied, 27,498 mapped untied, median de 0.0017.

In the control row there is no nearest-relative split, so "elsewhere" there means "mapped untied on a primary". The verdict uses the fractions only.

### 3. Selection caveats

- Every read of the 34-family and LRPAP1 sets was mapped on the combined primary assembly before this test, so UNMAPPED on `mat` or `pat` in the fate table means "mapped there, lost here" and was 0 in practice; reads unmapped on the primary assembly are the separate set of 959 (see above), and they enter no chain net (prereg S6 asked for unmapped reads in the net; this was not done).
- One animal and one tissue. The judged copies are 7 (mother) and 13 (father); no number here is a rate.
- 19 of the 20 judged copies have a relative within the allele cutoff on the reference, which limits what the verdict REFUTED says.
- The R4 read-level O2 readout (reads of a recovered copy before and after) was not computed.

### 4. Amendment 1 predictions

- **S1** (median de lower on the haplotype that has the copy, within the allele range): holds for 4 of 8 large loci with the mother as reference (GWFAM175_B0, 205_B0, 227_B0, p12) and 2 of 13 with the father as reference; the second clause (≤ 0.00958 on the present haplotype) holds for 20 of 21. Failures on the mother's side are copies with an identical or near-identical counterpart (equal medians) and GWFAM64_B11 (higher on `pat`). On the father's side many failures come from the label, not the biology: GWFAM70_B4's 97 reads sit at de 0.0013 on pat (the reference) and 0.0090 on the mat locus they are labelled to, so they belong to a copy the reference has.
- **S2** (the transcripts exist in both runs; kept as new and recovered only in the absent run): holds for GWFAM175_B0 (1 output on the locus in each run; kept new and recovered only in the `mat` run). For the other loci IsoCon has outputs on the locus in both runs (for example GWFAM208_B0 9 and 5) and none is kept as new, because they link to a reference locus within the cutoff.
- **S3** (the `pat` run applies the registered bars unchanged): reported above; R1 not testable, R2 and R3 pass, no candidate in a truth locus.

### 5. Arms

#### Reference = `mat`

| arm | R1 | R2 | R3 (bar ≤ 20%, registered rule) | R3 false-flag families, post hoc split | candidates / flagged | expressed copies recovered | flagged candidates by class |
|---|---|---|---|---|---|---|---|
| IsoCon | PASS (1 of 1) | PASS (0 of 13 flagged) | FAIL (6 of 24) | other-genome sequence only: GWFAM309, GWFAM318, GWFAM70, LRPAP1; matches nothing: GWFAM169, GWFAM382 | 21 / 9 | 1 of 16 | {'c_unmatched': 2, 'b_other': 6, 'a_recovered': 1} |
| in-house stage | FAIL (0 of 1) | PASS (0 of 13 flagged) | PASS (1 of 24) | other-genome sequence only: none; matches nothing: GWFAM382 | 20 / 1 | 0 of 16 | {'c_unmatched': 1} |

p12 line: IsoCon outputs with identity × coverage ≥ 0.999 on p12 (paternal) 1, on p14 (paternal) 1, on p14 (maternal) 1, on a p12 counterpart in `mat` 0. The p12 and p14 outputs differ by 4 bp, below the allele cutoff, so the link step treats p12 as an allele of p14 (designed-undetectable). The output counted on p14 (maternal) is a 2,562-bp transcript that is about 6% from both p14 copies (a different isoform) and is assigned to the maternal copy on 2 bp; no claim of three sequence types at one locus is made.

#### Reference = `pat` (held-out direction)

| arm | R1 | R2 | R3 (bar ≤ 20%, registered rule) | R3 false-flag families, post hoc split | candidates / flagged | expressed copies recovered | flagged candidates by class |
|---|---|---|---|---|---|---|---|
| IsoCon | NOT TESTABLE (0 of 0) | PASS (0 of 26 flagged) | PASS (2 of 21) | other-genome sequence only: GWFAM256; matches nothing: GWFAM175 | 22 / 7 | 0 of 30 | {'b_other': 4, 'c_unmatched': 3} |
| in-house stage | NOT TESTABLE (0 of 0) | PASS (0 of 26 flagged) | PASS (0 of 21) | other-genome sequence only: none; matches nothing: none | 16 / 1 | 0 of 30 | {'c_unmatched': 1} |

### 6. What did not go as predicted

1. The registered prediction "REFUTED at every LARGE locus" fails at GWFAM205_B0 (SUPPORTED) and GWFAM64_B11 (between bars) on `mat`, and at GWFAM390_B1 (SUPPORTED) and GWFAM382_B0 (PARTLY) on `pat`; and 19 of 20 judged copies are within the cutoff of a reference relative, so most of the REFUTED verdicts test nothing about hidden copies.
2. p12 is not a mother-absent locus (Amendment 2); the registered rule leaves the chr12 LRPAP1 difference unresolved, and the chain's transcripts point to c01.
3. The synteny-based truth rule misstates absence both ways: it calls copies with identical counterparts absent (GWFAM208_B0, 214_B0, 64_B17) and calls copies with no counterpart within the cutoff present (c01, GWFAM309/318/70 candidates), and `bonly` cannot list copies on chromosomes the primary assembly took from the reference haplotype.
4. S1 fails for 4 of 8 (mat) and 11 of 13 (pat) loci.
5. IsoCon R3 FAILS on `mat` by the registered rule (25%); the in-house stage fails R1.
6. The in-house arm is not the registered one: the binary contains the Amendment 15/15b fix (commit a13b817f), contrary to the prereg's "not applied" (Amendment 2); its acceptance re-run was not run.
7. Observation, not a result: GWFAM382_B0's 22 reads are 5 to 6% divergent from both haplotype assemblies, so they come from a copy present in neither genome, the class of the 2026-08-13 reference-absent candidates; this study did not test it.
8. The bar clause "with de ≈ divergence" is not in the verdict function.

### 7. Provenance

minimap2 2.30-r1287 (baseline BAM 2.31); reads aligned with `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M`; nearest-relative identities asm20 `-N 50 -p 0.1`; IsoCon 0.3.3 (`IsoCon pipeline --nr_cores 4`, defaults); `o3_candidates` binary built 2026-10-06 09:09 (sha256 prefix 1e01f514e387cadf), after commit a13b817f; delta 0.00958; the `mat` verdicts were written to `W/mat/VERDICTS.txt` at 2026-10-08 08:50:06 and the `pat` files were created after it, the `pat` run applying the same code and constants. Label cross-check against Amendment 10 (`express`, both haplotypes jointly): same order for GWFAM175_B0 304/281, GWFAM205_B0 82/77, GWFAM227_B0 21/13, but GWFAM208_B0 551/0, GWFAM214_B0 122/0, GWFAM64_B17 31/0, GWFAM64_B11 21/0 (those reads tie across haplotypes in Amendment 10 and are untied on `pat` here). Execution rulings and deferred minors are listed in section 9 and in Amendment 2.

### 8. Not covered

Between-individual differences (one animal); the testis library; copies absent from both haplotypes (GWFAM382_B0 above); the two `_pri` copies with a lift below 50% (GWFAM175:2, GWFAM491:1); which chr12 LRPAP1 locus is the extra one (c01 vs p12); a truth table that covers chromosomes the primary assembly took from the reference haplotype; unmapped reads in any chain net; the R4 read-level O2 readout; the Amendment 15 acceptance re-run of the in-house stage; any threshold tuning (none); independent replication.

### 9. Execution rulings and deferred minors (from the execution ledger, kept because the ledger is deleted)

#### Rulings
- Task 2: Ruling: tools/rlock.sh has no exec bit — invoke as 'bash tools/rlock.sh' everywhere (plan says tools/rlock.sh) — do not chmod a tracked file — cost if wrong: none
- Task 4: Ruling: p12 lifts as truth_lift class T? to mat chr12:27.43M (97.2% id, cov 1.0); chr12 cluster 4 pat loci vs 3 mat with many-to-one lifts (c01,c03 -> mat 34.29M), so the registered rule finds only chrY (sex) absent from mat for LRPAP1 — p12 kept as DESCRIPTIVE kind lrpap1_desc (no bar, not in R4 expressed), prereg Amendment 2 appended BEFORE labelling — cost if wrong: p12 would be reported as mat-absent when its mat counterpart may be an unrelated diverged copy; the user's LRPAP1 headline shrinks to 'chr12 has one more locus on pat'
- Task 4: Ruling: all 83 reads whose pat primary is on p12 are TIED on pat (p12 vs p14 pat, 99.84%) so the registered untied-placement label gives p12 0 reads — descriptive loci (lrpap1_desc) are labelled by the read's PRIMARY record on the truth haplotype (ties allowed), restricted to the LRPAP1 net; applies to no catalog/bar locus — cost if wrong: p12 reads may include reads that truly belong to p14 (the two are indistinguishable by score, which is itself reported)
- Task 6: Ruling: the mat FASTA order differs from the index/BAM-header order (name_map order-pairing assertion fired: JAQQLJ020000027.1 17,248 bp vs chr3_mat_hsa4_random_utig4-32 194,496 bp) but the length multisets are identical (225, one duplicated length) — pair by length (unique), duplicated lengths by relative order; refuse if the multisets differ — cost if wrong: two same-length scaffolds could swap names (no catalog locus or LRPAP1 copy lies on a scaffold; checked by the nonexistence of such hits in the loci tables)
- Task 6: Ruling: refabsent/copies.<hap>.paf targets are FASTA accessions (CM..), LRPAP1 PAFs use index names — panel chromosomes must be index names (BAM/FASTA/mmi): rename accessions through the same length-based name_map before building the panel (rename_hits) — cost if wrong: none (panel_to_copies fails loudly on an unknown name)
- Task 8: Ruling: o3_candidates on disk (built 2026-10-06 09:09, sha256 prefix 1e01f514e387cadf) post-dates commit a13b817f (2026-10-05 11:27 'Amendment 15/15b consensus fix — majority test, cs normalisation') so it INCLUDES the Amendment 15/15b fix, contradicting the prereg/plan text ('as shipped at b29afa55, Amendment 15 not applied'); use it as is (the fix repairs the documented A14 false-flag defect; a b29afa55 rebuild is a full cargo build and would test a known-broken stage), disclose in prereg Amendment 2 + artifact + write-up as 'o3_candidates with the Amendment 15/15b fix, acceptance (A13/A14/held-out re-run) NOT run' — cost if wrong: the in-house arm's false-flag behaviour is unvalidated rather than known-bad
- Task 12: Ruling: artifact page checked with ONE desktop-width light render (artifact-design rule: one look, then publish); dark theme and 400-px phone width were NOT rendered (palette validated in both modes with validate_palette.js; CSS built mobile-first with wrapping grids and overflow containers) — cost if wrong: a cosmetic issue on phone/dark that the user reports. Artifact: https://claude.ai/artifact/R1QoJa9GLfxnt51haKG6k8 (private). In-house caption states the Amendment 15/15b fix, not the plan's 'known defect'.
- Final: Ruling: the registered verdict function stays as implemented (fractions only) and all post hoc breakdowns are labelled post hoc — changing the verdict after seeing results would be tuning — cost if wrong: REFUTED labels stay generous for loci with identical counterparts, which the write-up now says
- Final: Ruling: reviewer 'Declined to judge' — (a) Amendment 2 ordering stands: it was appended before labels were computed (file edit precedes the labels run in this session; the same commit carries it with truth.py) ; (b) real divergence vs consensus error of the GWFAM309/318 candidates: unknown, stated as unknown; (c) in-house binary internals: accepted as is and disclosed; (d) artifact scan: stands

#### Deferred minors from the final review
- '8.7% away' for 175_B0 — 277/285 reads sit on the second asm20 hit (id .9108), not the registered nearest (.9126; 17 reads)
- 205_B0/390_B1 ties are between two equidistant within-cutoff copies (mechanism not stated)
- 'UNMAPPED = 0 by construction' is 0 observed (the chr12 reads show the converse is possible)
- GWFAM382_B0 de on pat 0.0565 (absorbed only) vs 0.0526 (all primaries) unexplained
- section 4 S1 'holds for exactly four' is post hoc
- label_reads family restriction not in prereg S5 (10 mat / 19 pat reads, none at a LARGE locus)
- score.recovered_loci does not filter kind (chrY/p12 hits count as recovered, never false); test pins it
- Review Focus tests 3 and 5 weak (verdict(None) unreachable; sex exclusion not distinguished from non-catalog)
- section 7 'first pat file 08:50:20' (lrpap1.genes.tsv 08:50:13; gate holds)

(Some minors were corrected in passing when the affected sentences were rewritten after the review: the unmapped-by-construction wording and the post hoc S1 wording.)


---

## PREREG_o3_maternal_reference_2026-10-08.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_maternal_reference_2026-10-08.md`)*

## Pre-registration — O3: what happens to the reads of copies the maternal haplotype lacks, and can a truth-free chain recover them?

**Written 2026-10-08, before any read is mapped or any number is looked at for this question.** DRAFT for user review; nothing below has been run.
Origin: advisor interest in the LRPAP1 page (UvFQf5QX): maternal and paternal haplotypes of one animal carry different copy numbers. Real-data version of
**Arm B** of `docs/archive/2026-09/PREREG_o3_reference_bias_2026-09-23.md` (its Arm A, simulation, refuted "missing copies hide as unmapped/tied reads" to 5% divergence).
Reuses, unchanged, the chain and bars of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` (Amendments 8, 9, 10). **No new thresholds.**

### 1. Questions
- **Q1 (fate).** Reads from a copy present on the paternal haplotype and absent from the maternal one: aligned to the maternal assembly, are they UNMAPPED, TIED
  (multi-mappers), or ABSORBED as untied primaries on another copy? Per locus, never pooled.
- **Q2 (O3 recovery).** Starting from the maternal alignment only (truth-free), does the net -> cluster -> flag chain recover those copies? Two clusterers,
  same inputs, same bars: IsoCon 0.3.3 and the in-house `o3_candidates` stage.

### 2. Substrate
KB3781 ("Jim", male), fibroblast Iso-Seq = the assembly's own cell line (`project_o3_matched_individual`). Reference = maternal assembly `mat`
(`GCA_028885495.2`, 225 sequences, `winloci_data/mGorGor1.mat.splice.mmi`); truth = paternal `pat` (`GCA_028885475.2`, 24 sequences). Testis (OR6737) is a
different animal: **excluded** (neither haplotype is its truth). Gorilla only; no human numbers.

### 3. Truth: MAT-ABSENT loci
A locus is MAT-ABSENT iff it exists on `pat` and has no counterpart on `mat`:
1. **Catalog loci**: `refabsent/bonly.tsv` rows with `hap = pat` (haplotype-only loci on the 9 chromosomes `_pri` took from `mat`; asm20, identity >= 0.90,
   coverage >= 0.80, overlapping no lifted copy of the family) — 35 loci / 34-family file `fams_bonly.txt`. Expressed = >= 3 reads best-placed there
   (`bonly_expressed.tsv`); LARGE = >= 20 reads. Known expressed: GWFAM175_B0 (281), GWFAM205_B0 (77), GWFAM227_B0 (13), GWFAM175_B1 (6).
2. **LRPAP1** (not in the catalog): the 11 loci of `lrpap1.copies.tsv` + `partial.copies.tsv`; mat-absent = those with no mat hit in `copies8.mat.paf` /
   `partial3.mat.paf`, read from the file and not from the artifact flags: **p12** (chr12 23.07 Mb, LOC134756368, 83 fibroblast reads) and the **chrY** copy.
3. **chrY copy = SEX CONTROL, not a copy-number locus** (male: Y is paternal only; same confound as the X-linked MAGE probes, `project_haplotype_cnv_proven`).
   Reported separately; excluded from every bar.
4. **Reverse direction** (held out; reference = `pat`, truth = `mat`): `bonly.tsv` rows with `hap = mat` (92 loci; expressed: GWFAM70_B4 24, GWFAM390_B1 20,
   GWFAM70_B1 16, GWFAM26_B2 13, GWFAM26_B3 9, ...). Primary direction is `mat` (the user's question); `pat` is the substrate held back.

SHARED-COPY CONTROL: reads whose best `pat` record falls on a pat locus that has a `mat` ortholog (lift in `copies_lift.tsv`, classes T2d/T2i).

### 4. Reads (selection bias disclosed)
- **R34**: `refabsent/scored.fa`, 32,219 net reads of the 34 families, already aligned to `mat` and `pat` with the baseline `@PG` command
  (`reads.{mat,pat}.bam`). Selection = a baseline (`_pri`) record on a copy of the family, so **every read was mapped on `_pri` before this test**.
- **R_LRP**: reads of the 11 LRPAP1 loci by the same rule (any primary/secondary record overlapping a locus body, <= 2,000 per locus, seed 1); aligned to
  `mat` and `pat` with the baseline `@PG` command (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, minimap2 2.30 vs 2.31, disclosed).
- **R_unm**: the unmapped primaries of the fibroblast BAM (959 reads, median 69 bp), aligned to `mat`/`pat`; descriptive only.
- Consequence: UNMAPPED on `mat` is observable only for reads that mapped on `_pri`; reads unmapped on `_pri` are covered by R_unm alone.

### 5. Read truth label and fate (Q1)
A read belongs to locus L iff its best `pat` record (alignment score, untied at 0.98 — the `express` rule of `refabsent_truth.py`) overlaps L.
Fate on `mat`, mutually exclusive, in this order (every constant is already registered):
1. **UNMAPPED**: no record, or primary query coverage < 0.80 (the registered hit-coverage rule; < 0.80 reported as PARTIAL, summed into this class).
2. **TIED**: primary MAPQ 0, or a secondary within 0.98 of the primary's alignment score (the registered tie rule).
3. **ABSORBED**: untied primary elsewhere; split into ON-NEAREST-PARALOG (the mat locus with the highest identity to L) and OTHER. Median `de` reported.
Per locus report: n reads, fate counts, nearest-mat-paralog identity (asm20 of L's sequence against `mat`), `de` of absorbed reads vs `de` of shared-copy reads.

**Bar (verbatim from the 09-23 registration, applied per LARGE locus, `pat`-copy reads on `mat`):**
>= 80% ABSORBED with `de` ~ divergence, <= 10% UNMAPPED, <= 10% TIED => advisor's claim **REFUTED at that locus**; 20-50% UNMAPPED+TIED => **PARTLY**;
>= 50% UNMAPPED or TIED => **SUPPORTED** at that locus.
**Prediction, written before looking:** REFUTED at every LARGE locus; GWFAM175_B0 ABSORBED at `de` ~ 0.066 (the `_pri` arm of Amendment 10 is this very
configuration: chr5 is mat-sourced); p12 ABSORBED on the chr14 ortholog of p14 (mRNA differs by 4 bp in 2,446; mat best hit NM 66 / 15.1 kb at MAPQ 60);
chrY copy ~ same fate as an autosomal absent copy of its nearest paralog. UNMAPPED > 10% is expected only where the nearest mat relative is far (see context panel).
**Context panel (no decision):** the 2026-08-14 whole-genome excision panel (`project_o3_excision_wholegenome`: ABSORBED 64.2%, ORPHANED 33.3%) re-plotted as
fate vs nearest-surviving-relative identity, so the natural cases (nearest relative >= 0.91) sit on the right of a spectrum whose left half is synthetic. Reuse of
existing tables only; if a per-family table cannot be recovered the panel is omitted, not recomputed.

### 6. Recovery chain (Q2) — reference = `mat`, truth-free until scoring
Net per family = reads with a record (primary or secondary) on any `mat` copy locus of the family (loci = `copies.mat.paf` hits at identity >= 0.90, coverage
>= 0.80; for LRPAP1 the 11 loci's `mat` hits), plus unmapped reads attributed as in A13; <= 1,000 reads per family, seed 1 (Amendments 9/10).
- **Arm I (IsoCon):** `IsoCon pipeline --nr_cores 4`, defaults, env `isocon` (as 10-01). Outputs -> flag (identity x coverage < 0.999 vs `mat`) -> link to a `mat`
  locus at whole-length d <= delta = 0.00958 (allele) -> merge (Amendment 8) -> candidates. **FLAG iff >= 2 transcripts** (floor registered in Amendment 9/10).
- **Arm H (in-house):** `o3_candidates` as shipped at main b29afa55, **opt-in, Amendment 15 fix NOT applied** (the A14 consensus defect is known and
  disclosed; this arm will probably reproduce its false flags — that is data, not a surprise). Same net, same link/merge/flag rules, same scorer.
- **Scoring (truth used only here):** a candidate RECOVERS locus L iff its best `pat` hit (identity x coverage >= 0.999) overlaps L (the Amendment 10 scorer,
  `refabsent_score.py`, with `mat` as baseline instead of `_pri`).

**Decisions (reused bars, applied per arm, mat-absent loci of Section 3):**
- **R1 (= D1)** every expressed BEYOND-delta LARGE locus is recovered by a FLAG. Smaller expressed loci reported, not judged.
- **R2 (= D2)** no expressed WITHIN-delta locus is flagged as a new copy — the delta rule's own check.
- **R3 (= D3)** families with no mat-absent expressed locus: those with a >= 2-transcript class b/c candidate <= 20%.
- **R4 (reported, not judged):** candidate-level sensitivity, precision and bipartite matching (`feedback_report_sens_prec_bipartite`); IsoCon outputs per
  recovered transcript; read-level O2 effect (reads of a recovered locus: primary on `mat` paralog vs on the candidate, `de` before/after, Amendment 10 O2 side).
- Neither arm is declared the winner; each bar is reported per arm.

**p12 is pre-registered as a DESIGNED-UNDETECTABLE case, not a failure.** Its mRNA is 4 bp (0.16%) from p14, below delta (0.96%), so the registered link rule
will classify it as an allele of the `mat` p14 locus. It is reported on its own line with three readouts: (a) do IsoCon / in-house separate p12 reads from p14
reads at the sequence level; (b) is the p12 transcript linked to p14 by delta; (c) exploratory **haplotype count**: distinct sequence types at the p12/p14
locus vs the 2 alleles one locus can carry (KB3781 has p14 on both haplotypes, so p12 + pat-p14 + mat-p14 = 3 types at one apparent locus). This rule is
**not** part of any decision here; if it looks useful it gets its own registration. Cluster-merging risk is stated up front: IsoCon merged TBC1D3 paralogs 6
edits apart in simulation (`RNA_ALLELE_ISOCON_2026-10-01.md`).

### 7. Outcomes that would change what we claim
- Q1: any LARGE locus with UNMAPPED+TIED >= 50% => the advisor's mechanism holds there; the O3 route would need an unmapped/tied-read recovery step (Arm B/C of 09-23).
- Q2: R1 fails for GWFAM175_B0 (recovered once before, 4 transcripts, id 0.9991) => the chain does not transfer from `_pri` to a haploid maternal reference.
- n is small by construction (4 expressed mat-absent catalog loci + p12): this is a demonstration per locus, not a rate. No pooled percentage will be quoted.

### 8. Not covered
Between-individual copy differences (one animal); the testis library; loci absent from both haplotypes; reads unmapped on `_pri` beyond R_unm; any threshold
tuning (none permitted); the Amendment 15 consensus fix (separate registered work).

### 9. Compute plan (CRASH RULE: foreground, serial, one heavy job at a time; `/mnt/linuxdisk` has ~25 GB free, so no new index is written)
(1) extract R_LRP + R_unm; (2) map to `mat` then `pat` (existing .mmi, ~14 GB RSS each); (3) fate tables; (4) nets; (5) IsoCon per family (resumable
`iso_batch.sh`); (6) in-house `o3_candidates` per family; (7) link/merge/score; (8) artifact. Register rows to be drafted after results, not before.

---

### Amendment 1 (2026-10-08, written before any read is mapped or any number is looked at) — paternal-reference run and side-by-side

DRAFT for user review. Requested by the user after approving the sections above: show side by side what the same copies look like where they
are PRESENT (the paternal genome) and where they are MISSING (the maternal genome), and how the recovery chain's work in the maternal run compares
with the paternal run, where the same transcripts are simply in the reference. Changes S3 item 4 and S6; everything else stands.

**A1.1 Symmetric runs.** S3-S6 are run twice with one switch (`O3_REF`): reference `mat` (primary, exactly as registered above) and reference `pat`
(truth = `mat`). Same read sets (R34, R_LRP, R_unm), same indexes, same 34 families + LRPAP1, same constants, same bars; nothing is tuned between runs.
**Order:** the `mat` run is completed and its verdicts written to `W/mat/VERDICTS.txt` before the `pat` run starts (a rule validated where it was
developed is not validated: the `pat` run is the substrate held back).
- Truth for `O3_REF=pat`: `bonly.tsv` rows with `hap = mat` (92 loci; known expressed: GWFAM70_B4 24, GWFAM390_B1 20, GWFAM70_B1 16, GWFAM26_B2 13,
  GWFAM26_B3 9, ...), plus LRPAP1 loci present on `mat` and absent from `pat` (expected none: the PAFs decide, as for `mat`).
- Supersedes S3 item 4 ("held out, not run"): the reverse direction is run once, after the `mat` run, and reported per locus.
- Not covered: the two `_pri` copies on mat-sourced chromosomes with lift < 50% that `truth_classes` never refined (GWFAM175:2, GWFAM491:1) are
  not in the truth set of either run.

**A1.2 Side-by-side readouts (descriptive; no new bar).** For every absent locus with >= 3 reads, per run:
- (a) *paired reads*: for each read, `de` and MAPQ of its primary on the reference haplotype (copy missing) and on the other haplotype (copy present);
- (b) *divergence pile*, LARGE loci only: per-bin primary-alignment coverage and mismatch counts (60 bins, positions as fractions of the interval) along
  the nearest reference paralog (where the reads land) and along the locus itself;
- (c) *the chain's view of the same copy in the two runs*: IsoCon outputs whose best hit (identity x coverage >= 0.999) overlaps the locus, counted
  against the truth haplotype in run A and against the reference in run B; how many the link step kept as new copies in each run; the flagged
  candidates recovering the locus (run A); the in-house stage candidates whose nearest locus overlaps it (run B: `d`, flagged).
(a) is partly circular on the `other` haplotype by construction (a read is labelled to a locus by its untied best placement there), so the
`other`-side fate is not reported as a result; only `de`, MAPQ and the count of reads dropped as ambiguous are.

**A1.3 Predictions, before looking.**
- **S1** For every LARGE locus the median `de` on the haplotype that has the copy is smaller than on the haplotype that lacks it; on the present
  haplotype it lies within the registered within-species range (<= delta = 0.00958).
- **S2** The IsoCon outputs of an expressed absent copy exist in both runs (they come from the reads, not from the reference): in run B they match the
  reference at >= 0.999 and none is kept as new (the flag step skips such outputs by construction — shown for consistency, not as a test); in run A they
  are kept as new and, where R1 passes, recovered.
- **S3** In the `pat` run the registered bars R1-R3 are evaluated unchanged for the mother-only copies. No outcome is predicted beyond the `mat` run's;
  it is reported even if it fails.

**A1.4 Visual contract.** The artifact shows both directions behind a toggle; per LARGE locus a left/right pair (reference where the copy is missing |
other haplotype where it is present) with the paired-read plot, the divergence pile, the fate bars and the chain's view.

**A1.5 Cost.** About doubles the compute of S6 (second IsoCon pass, second set of in-house batches, second set of contig alignments). No new index.

---

### Amendment 2 (2026-10-08, written during Task 4 after the truth loci were built and BEFORE any read was labelled or mapped to a locus) — p12 is not a mat-absent locus by the registered rule

Finding (from the PAFs the prereg names as the authority, not from the artifact flags): `truth_lift` classes p12 `T?` (lift 0.914, 112 mismatches) and lifts it to
`CM054594.2:27,426,744-27,441,619`, where `mat` carries a locus matching p12's body at identity 0.972 / coverage 1.0 (asm20 -N 50 -p 0.5). The chr12 LRPAP1 cluster
has 4 loci on `pat` (22.55, 23.07, 24.79, 30.2 Mb) and 3 on `mat` (25.87, 27.43, 34.29 Mb), but the lifts are many-to-one (c01 and c03 both lift to mat 34.29 Mb, c01 with
222 mismatches), so the registered synteny rule cannot say WHICH `pat` locus has no `mat` counterpart. Consequence for S3.2 / S5: under the registered rule the only
LRPAP1 locus absent from `mat` is the chrY copy (sex control); p12 is NOT a mat-absent locus.
Ruling (no new threshold): p12 is kept as a DESCRIPTIVE locus, kind `lrpap1_desc` = an LRPAP1 locus on a chromosome `_pri` took from the truth haplotype whose lift
class is `T?`. It gets a fate row, side-by-side panels and the p12 chain line; it carries NO bar verdict, is not counted in R1-R4's expressed set, and is labelled
in the artifact as "diverged counterpart on the reference". The S5 prediction for p12 ("absorbed on the chr14 ortholog of p14") stands as a prediction about where
its reads go; it is no longer a test of the mat-absent claim. The copy-number statement for LRPAP1 on chr12 is "4 loci on pat, 3 on mat; which one is extra is not
resolved by the registered rule".
Amendment 2, labelling note (same day, before any fate was computed): on `pat` all 83 reads whose primary record lies on p12 are TIED with p14 (score ratio within 0.98; the two
copies differ by 4 bp in 2,446), so the registered untied-placement label gives p12 zero reads. Descriptive loci (`lrpap1_desc`) are therefore labelled by the read's
PRIMARY record on the truth haplotype, ties allowed, restricted to the LRPAP1 net. No catalog locus and no bar uses this rule; the tie itself is reported as a result.
Amendment 2, arm H identity (same day, before the in-house stage was run): S6 says the in-house arm is `o3_candidates` "as shipped at main b29afa55, Amendment 15 not applied". The
repository history shows the Amendment 15/15b consensus fix (majority test, cs normalisation; commit a13b817f, 2026-10-05) was implemented after b29afa55, and the binary on disk
(`rustle_target/release/o3_candidates`, built 2026-10-06 09:09, before the 12:01 source consolidation) post-dates it. Arm H is therefore **`o3_candidates` WITH the Amendment 15/15b
fix**, whose acceptance re-run (A13/A14/held-out) was not run. Its flags are reported as such; the "known false-flag defect" caption in the plan no longer applies and is replaced by
"Amendment 15/15b fix present, acceptance not run".


---

## PREREG_unmapped_rescue_2026-10-08.md

*(original: `docs/archive/2026-10/merged/PREREG_unmapped_rescue_2026-10-08.md`)*

## Pre-registration — rescuing reads that became unmapped: a pool-first, reference-free net

**Written 2026-10-08, before any clustering or attribution result exists for this question.** Goal (user, 2026-10-08): design a better net or clustering that saves unmapped reads, either
those that are unmapped on one parental haplotype or those left unmapped after a copy is erased from the genome and the reads are realigned.
Standing context: `docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md` (rows 1217, 1221, 1223): after the A13 net (attribution of EACH unmapped read by alignment to the family's reads and copies, read coverage >= 0.5, de <= 0.20)
23 of 53 deleted copies still had no read in their family's net and 7,184 deleted-copy reads were out of reach; 31-mers placed 0 of 5,312 unmapped reads. `docs/NEGATIVE_RESULTS_REGISTER.md` row 1065: a pile of unmapped reads is emitted as a FLAGGED CANDIDATE, never a family;
IG/TR and run-exclusivity screens are compulsory for anything called a new copy. This study rescues READS into a net; it does not define families.

### 1. The idea (frozen)

The A13 net judges every unmapped read alone against the family. Reads of one deleted copy are near-identical to EACH OTHER (one transcript, HiFi errors) even when each is 15-30% from every surviving relative. So:
1. **Pool first.** All-vs-all overlap of the unmapped pool alone (no reference). A read-read edge exists iff the overlap's gap-compressed divergence `de` <= delta = 0.00958 (the registered allele cutoff) and the aligned block covers >= 50% of the shorter read (the registered merge rule of Amendment 8).
   Clusters = connected components; a cluster needs >= 3 reads (the stage's `--min-cluster 3`).
2. **Consensus.** One abPOA consensus per cluster (<= 100 reads, longest first, seed 1).
3. **Attribute the cluster, not the read.** Translated search of the consensus (mmseqs2, `--search-type 2`, both strands, e-value 1e-5 as the external BLAST convention) against the surviving copies of the family catalog (915 copies minus the erased ones, DNA of the unmasked `_pri` intervals).
   The cluster goes to the family of its best hit iff that family's best bit score is strictly above every other family's (no tie); otherwise the cluster ABSTAINS.
4. **Rescue.** Every read of an attributed cluster joins that family's net.
Comparators: B0 = A13's single-read rule re-run on the same pool (minimap2 map-ont to the surviving copies, read coverage >= 0.5, de <= 0.20); B1 = single-read translated attribution (the same search as step 3, one read at a time, same tie rule). B1 isolates what pooling and consensus add.

### 2. Beds and what is held back

| bed | substrate | pool | truth | role |
|---|---|---|---|---|
| **A** (dev) | 2026-08-14 whole-genome excision, 162 two-copy families, `o3_excise/masked.bam` (347,798 reads, masked genome) | the 61,512 unmapped records + the 959 unmapped primaries of the fibroblast BAM (background) | a read is D of family F iff its baseline primary (`baseline.tsv`) lies in the erased interval of F | all design checks and the sensitivity sweep are on A |
| **H** (held out) | A13's 53 held-out multi-copy families, `rna_allele/linktest/R.bam` (masked), `labels.tsv` role D | its 5,313 unmapped records + the same 959 background reads | `labels.tsv` (copy of origin) | **not looked at until the rules in section 3 are frozen by a commit**; run ONCE |
| **M** (real) | the 959 unmapped primaries of the KB3781 fibroblast BAM against the maternal haplotype (`o3_mat/map/R_unm.mat.bam`: 132 map, 125 of them in one 71-kb stretch of mat chr12 95.96-96.03 Mb) | those 959 reads | the maternal alignment | held out like H: run ONCE, no tuning |

### 3. Rules frozen before the first run on any bed

Decided on bed A only, written here, committed, then H and M are run once with no change. Constants in section 1 are the only ones; the sweep on A reports delta/2, delta and 2 delta for the edge rule and e-values 1e-3 / 1e-5 / 1e-10, and the registered setting is delta and 1e-5.

### 4. Metrics (per bed, per method; read level unless stated; never pooled across beds)

- N_D = unmapped D reads; coverage = share of N_D in a cluster of >= 3; purity = share of clustered D reads whose cluster's majority truth family equals their own; fragmentation = clusters per deleted copy.
- Attribution: attributed / abstained clusters; cluster accuracy = attributed clusters whose majority family is the attributed family.
- **Rescued-correct** = D reads of an attributed cluster attributed to their own family; **wrong joins** = reads of an attributed cluster attributed to another family, and background reads in any attributed cluster.
- Deleted copies with >= 1 rescued-correct read (copies reached) against the copies that have unmapped reads.
- Candidate fidelity: the consensus of each attributed cluster aligned to the erased copy's transcript sequence (asm20): identity x coverage >= 0.999 counts as the copy reproduced from the unmapped reads alone.
- Bed M: whether the chr12 reads form one cluster, the consensus' identity to the maternal region, and what the 834 other reads do.

### 5. Bars (committed now)

- **U1 (the pool is coherent)** on H: coverage >= 0.90 and purity >= 0.99.
- **U2 (the net is better)** on H: rescued-correct >= 2,656 (half of the 5,312 unmapped D reads; A13 rescued 557) with wrong joins <= 5% of rescued reads (A13's registered false-move bound is 5%; it measured 0.77%), AND rescued-correct of the pooled method > B1 > B0 on the same pool.
- **U3 (real)** on M: >= 90% of the 125 chr12 reads in one cluster, consensus identity to the maternal region >= 0.99, and no background read in that cluster.
- Reported, not judged: copies reached, candidate fidelity, behaviour of the sweep, runtime.

### 6. Predictions, before looking

U1 holds (components at delta over reads of one transcript are nearly pure; the risk is fragmentation, not mixing). U2 is uncertain: attribution decides it, and translated search should reach coding families but not pseudogenes, lncRNA or very divergent copies;
I expect the pooled method to beat B1 and B1 to beat B0, and rescued-correct between 20% and 60% of N_D, so U2's 50% bar may FAIL. U3 holds for the cluster (one expressed sequence, 125 reads, 3.6 kb) with no family to attribute it to.

### 7. What would change the claim

If clustering holds (U1) but attribution fails (U2), the supported statement is "the unmapped reads of an erased copy are recoverable as coherent candidate transcripts (flagged, not family-attributed)", which is the register's own allowed class; the family label then needs the annotation-free evidence of the other layers.
If U1 fails on H, the pooled net is refuted as designed.

### 8. Not covered

Reads that map poorly (absorbed on a paralog) rather than unmap; non-excision absences (a copy real in another individual); the full 23 GB library as one pool (cost); any change to the shipped binary; IG/TR and run-exclusivity screens beyond reporting them for flagged clusters.

### 9. Execution

Python on top of minimap2 (`-x ava-pb`), pyabpoa, mmseqs2 and networkx, in `bench/unmapped_rescue/` with unit tests for the graph, the consensus selection, the tie rule and the scoring; heavy jobs foreground and serial under `tools/rlock.sh heavy` (`bash tools/rlock.sh`); work dir `/mnt/linuxdisk/tmp/o3_rescue/`.

---

### Amendment 1 (2026-10-08, after the dev-bed A results below, before bed H or any bed-M-specific step is run)

**Seen on bed A before this amendment (dev):** the seed-round clustering at the registered delta puts 57,664 of the 59,791 deleted-copy reads (96.4%) in clusters of >= 3 at purity 1.0000 (221 clusters, 1 to 7 per family, median 2).
Mmseqs translated attribution of the 221 consensus sequences (section 1 step 3 as written) rescued 77.8% of the reads but 15% of the joined reads were wrong joins and cluster accuracy was 0.59; the wrong clusters were near-ties between paralogous families (median best / runner-up bit-score ratio 1.03).
Replacing the score by the NUCLEOTIDE evidence of the whole consensus changed this: with discontiguous megablast (`blastn -task dc-megablast`, e-value 1e-5) the family score = the number of consensus bases covered by the union of the HSPs of the family's best target, and the dev curve is flat:
strictly-above rule 87.8% rescued / 7.1% wrong joins; relative margin >= 0.10 87.8% / 0.1%; 0.20 87.8% / 0.1%; 0.50 87.7% / 0.0%; an absolute margin of 10 bp 87.8% / 0.6%; a floor on the covered fraction of the consensus costs reads without helping (0.20: 68.3% / 0.7%).

**Frozen for beds H and M (replaces section 1 step 3; everything else stands):** attribute a cluster iff its best family's covered-base score is >= 1.10 x the runner-up family's score (a family with no hit scores 0, so a single-family hit is attributed); otherwise abstain. The 1.10 comes from bed A where 0.10 to 0.50 give the same answer; H and M report the sweep 1.00 / 1.10 / 1.50 and the translated mmseqs score as a comparator.
Comparators B0 (A13's single-read rule) and B1 (single-read dc-megablast, same margin rule) are run on a seeded sample of 1,000 unmapped deleted-copy reads per bed (cost); the pooled method's numbers are on all reads.

**Bed M is not fully held out.** The 959 background reads are part of bed A's pool, so the clustering outcome for the 125 chr12 reads was visible in the dev run. U3's clustering part is therefore reported as dev-visible; its consensus-fidelity part (identity of the consensus to the maternal region) and the attribution outcome are first evaluated after this amendment.

**Additional readouts (reported, not judged):** deleted copies reached; per-family cluster fragmentation; the clusters that abstain (their reads stay unrescued, with size); whether the consensus reproduces the erased copy (identity x coverage >= 0.999 against the unmasked genome at the erased interval).

---

### Amendment 2 (2026-10-08, after the frozen rule had been run once on bed H; written BEFORE any rescue number of the new score exists)

**What was seen on H and why this is exploratory.** Asked whether a family profile would resolve the attribution failures, I (a) scored families by the bases covered by ANY surviving member (a profile proxy): identical to the frozen rule on H (4,255 rescued, 3 wrong, 20/32 copies) and on A (52,493 vs 52,496); and (b) searched each of the 13 abstained H clusters that has a truth family against ONLY its own family's survivors with permissive settings
(blastn word 7, e-value 10; tblastx): 11 of 13 show homology far above the random-family background (for example 475 bits, E = 5e-135, for the 200-read GWFAM158 cluster; random families <= ~150 bits), two do not (GWFAM244: 31 bits; one GWFAM348 cluster: 35 bits).
So the signal exists and the frozen score ranks it badly. This diagnosis used the truth of bed H, so **H is no longer held out for anything derived from it**; the rule below is evaluated first on dev bed A (whose abstained clusters were not inspected), then on H as an exploratory arm and reported as such.

**New score (parameter-free; the margin rule is unchanged).** A family covers a consensus base if any of its copies has a dc-megablast HSP over it; a base covered by n families is credited 1/n to each. The cluster goes to the best family iff its score >= 1.10 x the runner-up (the frozen margin), else it abstains.
Hypothesis: repeat-derived hits give unrelated families near-equal cover, which the frozen score cannot separate from the true family's exon homology; crediting shared bases less should raise the true family's rank.
**Decision (before looking):** the arm is judged by wrong joins <= 5% of rescued reads AND rescued-correct above the frozen rule on A (52,496) and on H (4,255). If it only helps on H it is reported as unsupported (tuned on the held-out bed). The profile-HMM arm (nhmmer on a per-family MSA) is NOT run unless this arm shows the signal is rank-limited rather than absent.

**Amendment 2, result on dev bed A (the arm as registered):** 52,516 rescued (87.8%, +20 reads) but 4,041 wrong joins (7.15% of joined, cluster accuracy 0.667): it FAILS the <= 5% condition. Cause (55 clusters changed, 4 newly right, 51 newly wrong): crediting shared bases 1/n shrinks every score, so noise-level top hits such as 13 against 7 credited bases pass the 1.10 margin. The margin is relative; the arm has no evidence floor.

---

### Amendment 3 (2026-10-08, before any floor number exists)

Add an evidence floor to the Amendment 2 score: the cluster is attributed iff its best family's specific cover score is >= f x the consensus length AND >= 1.10 x the runner-up. **f is selected on dev bed A only:** the smallest f in {0.02, 0.05, 0.10, 0.20} with wrong joins <= 1% of joined reads; if none qualifies the arm is dropped. The selected f is then applied once to H (exploratory: H is no longer held out, see Amendment 2) with the same conditions as the frozen rule (wrong joins <= 5% of rescued; rescued-correct above 4,255).
If the arm beats the frozen rule on both A (52,496) and H (4,255) the claim is "a repeat-robust cover score with an evidence floor attributes more clusters"; if it helps only on H it is unsupported. The per-family HMM profile arm stays unrun unless this arm leaves a rank-limited remainder on A.

**Amendment 2 and 3 results, corrected after the cap fix (see Amendment 4).** With BLAST targets uncapped the arms read: Amendment 2 (parameter-free specific cover, margin 1.10): dev A 55,934 rescued (frozen 55,901), 610 wrong (1.08% of joined); held-out H 4,325 rescued (frozen 4,255), 457 wrong (9.56%).
It passes the <= 5% condition on A and FAILS it on H. Amendment 3 (selection rule: smallest f with wrong <= 1% on A): f = 0.02 (A 270 wrong, 0.48%); on H it gives 4,325 rescued and 271 wrong (5.90% > 5%): FAILS on H. The earlier figures in the Amendment 2 and 3 text (7.15% wrong on A, the 51 newly wrong clusters) were produced with the capped BLAST output.
Post hoc, not selected by the rule: f = 0.05 gives A 55,852 / 14 wrong and H 4,321 / 27 wrong (0.62%), better than frozen on H and worse on A; it is reported, not claimed.

---

### Amendment 4 (2026-10-08, disclosure after the independent review; no rule is changed)

1. **BLAST cap.** The dc-megablast runs used `-max_target_seqs 500`, below the catalog size (753 / 862 / 915); BLAST truncates in a non-score order, so a family's own copy can be dropped (87 of 188 dev consensus sequences and 10 of 39 held-out ones reached the cap). Re-run with 5,000 and a guard (`attribute.capped_queries`):
   bed H is unchanged; dev bed A changes (frozen rule 52,496 -> 55,901 rescued, 87.8% -> 93.5%; strict rule 7.1% -> 1.1% wrong joins, so the margin mattered less on A than Amendment 1 stated). The capped outputs are kept as `*.cap500.*`.
2. **Deviations never amended:** the clustering is seed rounds (`map-hifi` against 2,000 seeds, <= 6 rounds) instead of the all-vs-all of section 1, because the all-vs-all did not finish; the comparator sample is 500 reads per bed, not 1,000; the code implementing the frozen rule was committed after bed H was run (the dev numbers of Amendment 1 reproduce with it); the delta sweep and bed A's rescue file were produced after bed H.
3. **Bed M target set.** The M verdict in the first write-up used bed A's target set (162 families erased). In M's world nothing is erased, so the target set is the full 915-copy catalog.
4. **B0** was specified in the body as the single-read comparator against the surviving genomic copies; it is not A13's rule (A13 also aligned to the family's net reads), so B1 > B0 is vacuous.
5. **Fidelity metric.** The registered metric is identity x coverage >= 0.999; the first write-up reported identity alone.

---

### Amendment 5 (2026-10-08, before any profile number exists): the per-family HMM profile arm (exploratory; bed H is not held out)

**Why now.** After the cap fix, the abstentions that a better model could win are the rank-limited ones: on H the true family ranks below the top two for 6 clusters (228 reads), second for 2 (94) and has no hit while another family does for 3 (89): 411 reads at most (7.7% of N_D); on the two-copy dev bed a profile is a single sequence, so only H can test it.
**Arm.** Universe = the 53 H families (the only catalog families whose surviving copies have mapped reads in the masked run). For each survivor, the transcript consensus = abPOA of its <= 100 longest primary-aligned reads; per family an MSA (MAFFT) of its survivors' consensus sequences and `hmmbuild --dna` (HMMER 3.4 from the `hmmer` mamba env); `nhmmer` of each family HMM against the cluster consensus sequences (e-value <= 1e-3).
Score = the consensus bases covered by the family's nhmmer hits (union of the alignment spans); attribute iff best >= 1.10 x runner-up (the frozen margin). Variant reported, not judged: the best bit score with the same margin.
**Fair comparison.** The frozen rule re-evaluated in the SAME universe (dc-megablast HSPs restricted to the survivors of the 53 H families). The HMM arm is a gain only if it beats that restricted baseline in rescued-correct with wrong joins <= 5% of rescued reads; if a family has a single survivor its profile is that sequence and no gain is expected, so the result is also reported for families with >= 3 survivors separately.
**Not a held-out claim.** The abstained clusters were inspected on H, so a gain here is hypothesis-generating; it would need a new multi-copy erasure bed (erase one copy of other multi-copy catalog families) to be confirmed.

**Amendment 5 result (H; exploratory).** HMM profile: 4,522 rescued, 125 wrong (2.69%), 22 of 24 copies; restricted dc-megablast baseline in the same universe: 4,279 / 59 (1.36%) / 21 of 24; bit-score variant 4,479 / 183 wrong. The condition (more rescued with wrong joins <= 5%) is met: +243 reads, from two clusters (families with 9 and 5 survivors).
Costs: the GWFAM163 profile attracts five wrong small clusters; one cluster lost. Families with >= 3 survivors: 496 rescued / 14.3% wrong vs 226 / 20.7%. Not a held-out result.

---

### Amendment 6 (2026-10-08, after a stop-hook review of the claim; written before any hybrid number exists)

**The claim, stated against the existing approach.** The existing net on bed H is A13's (register rows 1221 and 1223): 557 of the 5,312 unmapped deleted-copy reads joined, 99.6% right. The pool-first net with the frozen rule rescues 4,255 (80.1%) of the same reads, 1 wrong cluster of 27, with the same truth. B1 (single-read dc-megablast with the same score and margin) is a variant of the NEW design, not the existing approach;
it shows the gain over A13 comes from the sensitive nucleotide cover score with a margin, and that pooling adds cost (about 6x less wall-clock) and consensus transcripts but no extra rescued reads (406 vs 410 of 500).
**Hybrid (parameter-free).** Pool first (clusters, consensus, cluster-level attribution as frozen); then every read NOT rescued by an attributed cluster (reads of abstaining clusters, reads in no cluster, background reads) is attributed alone with the same cover score and margin 1.10 (single-read dc-megablast against the same targets, cap 5000).
**Reading, read level, per bed (A dev first, then H):** rescued-correct = deleted-copy reads attributed to their own family; wrong = deleted-copy reads attributed to another family plus any background read attributed to a family; cost = the wall-clock of the extra per-read alignments.
**Decision before looking:** the hybrid supports "pooling does not have to cost rescue" only if on BOTH beds its rescued-correct >= the pooled frozen rule's and >= B1's estimate on the 500-read samples scaled to the bed, with wrong joins <= 5% of joined reads, and its extra alignment work is <= 25% of per-read alignment of the whole pool.
If it does not meet this, the pool-first net is reported as better than A13 and cheaper than per-read attribution but not as a net that rescues more than per-read attribution.

---

### Amendment 7 (2026-10-08, at the user's request, before any number): the reconstructed loci as extra reference (augmented-reference realignment) on held-out bed H

**Idea (user).** Use the loci reconstructed from the reads (the cluster consensus sequences) as additional reference and realign the reads with minimap2, so each read lands on the reconstructed locus (and the consensus can be placed among the family's nodes). No family attribution is needed for this step.
**Design.** Reference = the masked genome (the existing baseline alignments `linktest/R.bam`, `splice:hq -uf --eqx -Y -N 50 -p 0.1`) PLUS the consensus sequences. Instead of a new 13 GB index, the reads are aligned with the SAME minimap2 command to the consensus sequences alone; per read the primary becomes the consensus record iff its alignment score `AS` is strictly higher than the genome primary's (ties stay on the genome; an unmapped read takes any consensus record at query coverage >= 0.80). This is the primary a concatenated reference would give; MAPQ is not reproduced.
**Non-circular test (primary).** Half of the unmapped deleted-copy reads (seeded split, seed 1) build the clusters and consensus sequences (background reads stay in the building half); the OTHER half, which never entered clustering, and the reads that were absorbed on a paralog in the masked run (deleted-copy reads with a genome primary), and the survivors' reads (S) are realigned.
A full-pool version (consensus built from all unmapped reads, mapping the same reads back) is reported as a self-consistency check, not as a result.
**Readings.** (a) held-out unmapped D reads: share that get a consensus primary at identity >= 0.98 whose cluster's majority family is the read's own family; (b) absorbed D reads: share that move to an own-family consensus, and median `de` before and after; (c) FALSE MOVES: S reads whose primary becomes any consensus, and background reads mapping to a family-cluster consensus.
**Bars (before looking):** (c) false moves of S reads <= 5% (A13's bound; A13 measured 0.77%) and background reads <= 5%; (a) >= 80% of held-out unmapped D reads that belong to a family with a consensus in the building half (the clustering coverage on H was 95.7%, so a miss is mostly a copy with no cluster in the half); (b) reported, not judged.
**Prediction:** (a) holds near 90%; (b) moves a minority of absorbed reads (those at the distant paralogs; reads on a 1% paralog stay unless the consensus is closer); (c) holds below 1%.

### Amendment 8 (2026-10-08, at the user's request, before any LRPAP1 number): LRPAP1 and its gorilla copies as the worked example

**Substrate.** KB3781 fibroblast reads of the 11 LRPAP1 loci (2,875 net reads, `o3_mat/map/R_LRP.{mat,pat}.bam`); reference = the maternal haplotype; truth = the paternal one. The LRPAP1 copies are 96 to 99% identical to each other, the regime where a missing copy hides as an absorbed pile (small `de`), not as unmapped reads.
**Pipeline (the same truth-free steps).** (1) the net = LRPAP1 reads whose primary on `mat` has `de` > 0.00958 (the registered allele cutoff) or is unmapped; (2) seed-round clustering of the net (same edge rule), consensus per cluster; (3) align each consensus to `mat` and `pat` (identity x coverage; which locus it lands on, in particular pat chr12 22.55 Mb, the copy IsoCon's transcripts matched); (4) augmented reference as in Amendment 7 (mat primaries from `R_LRP.mat.bam` plus the consensus sequences): which reads move, their `de` before and after, false moves of reads of other copies;
(5) the consensus placed among the family's nodes: its identity to each LRPAP1 copy transcript on `mat` (the copies' RefSeq/Liftoff intervals), nearest node first.
**Readings and bar.** Report each step. Bar for the example to count as a demonstration: a cluster whose consensus matches a pat-only locus at identity x coverage >= 0.999 and a mat locus at < 0.999 (the O3 flag), AND false moves <= 5% of the reads that were on the other copies.

#### Amendment 6 result (bed H, 2026-10-08; bed A not run, see below)

| H, 5,313 deleted-copy reads | rescued correct | wrong | unattributed | wrong / joined |
|---|---|---|---|---|
| pooled frozen rule | 4,255 (80.1%) | 3 | 1,055 | 0.07% |
| hybrid (residual reads one by one) | 4,421 (83.2%) | 71 | 864 | 1.58% |
| B1 estimate (410/500 scaled) | 4,357 (82.0%) | | | |

Residual reads: 2,014 of 6,272 (32.1%); extra BLAST 506 s (about 0.25 s per read, so per-read attribution of the whole pool is about 1,580 s).
**Decision rule:** rescued-correct >= pooled (4,421 >= 4,255, met), >= B1 scaled (4,421 >= 4,357, met), wrong <= 5% of joined (1.58%, met), extra alignment work <= 25% of the whole pool (32.1%, **NOT met**).
**Verdict: the hybrid does not meet the registered rule on H.** The pool-first net is reported as better than A13 and cheaper than per-read attribution, not as a net that rescues more than per-read attribution. The hybrid's +166 reads over the pooled rule and +64 over B1 are a descriptive observation, not a claim.
Bed A was not run: the rule needs BOTH beds, so the H failure fixes the verdict and a dev run cannot change it (disclosed; an hour of BLAST saved).

#### Amendment 7 result (bed H, 2026-10-08, `bench/unmapped_rescue/run_augment.py`)

Building half: 2,656 of the 5,313 unmapped deleted-copy reads (seed 1) + the 959 background reads, clustered with the frozen rule (coverage 95.1%, purity 1.0000, 44 clusters, 0 background reads in a family cluster), abPOA consensus per cluster. The other 2,657 unmapped D reads never entered clustering.
Reads are aligned to the consensus sequences with the baseline command; a read moves iff the consensus alignment score is strictly higher than its genome primary's (unmapped reads: identity >= 0.98, query cover >= 0.80).

| class (non-circular half run) | n | moved | to own family | median de before -> after |
|---|---|---|---|---|
| held-out unmapped D | 2,657 | 2,505 (94.3%; 94.6% of the 2,648 whose family has a consensus) | 2,505 | unmapped -> 0.0011 |
| absorbed D (genome primary on a paralog) | 11,973 | 2,250 (18.8%; 82.4% of the 2,728 whose family has a consensus) | 2,247 | 0.1237 -> 0.0031 |
| S (survivor copies' reads) | 41,727 | 0 (0.00%) | | |
| background | 959 | 142 (all onto background-only clusters built from them), 0 onto a family cluster | | |

Bars: (a) >= 80% of held-out unmapped D with a family consensus: 94.6%, **met**. (c) false moves of S <= 5%: 0.00%, **met**; background onto a family cluster <= 5%: 0.0%, **met**. (b) reported: absorbed reads move only when the consensus is closer (median de 12.4% -> 0.3%), 99.9% of the moves onto the own family.
Full-pool self-consistency (47 consensus from all unmapped D, same reads mapped back): unmapped D 4,973 of 5,313 (93.6%), absorbed D 2,248 (2,243 own family), S 2 of 41,727 (0.005%, both own-family consensus), background 0 onto a family cluster. Not a result; shown to confirm the half run loses little.
Prediction check: (a) predicted near 90%, observed 94.6%; (b) predicted a minority, observed 18.8% (82% in families with a consensus); (c) predicted below 1%, observed 0.00%.
Caveats: the score compares a spliced genome alignment with an unspliced transcript alignment; S moving 0 of 41,727 is the control that the comparison does not favor the consensus by construction. A concatenated reference was not built (no 13 GB index); the primary is the one it would give, MAPQ not reproduced.

#### Amendment 8 result (LRPAP1, KB3781 mother as reference, 2026-10-08, `bench/unmapped_rescue/run_lrpap1.py`)

**Net.** 381 of the 2,875 LRPAP1 reads have maternal divergence `de` > 0.00958 (none unmapped; 138 on the paternal haplotype). Clustering with the frozen edge rule (all net reads were seeds, so this is one all-vs-all round): 4 clusters, 311 reads, purity 1.000 by the paternal primary's locus. Cluster sizes 283 (c01), 21 (c00), 4 (c01), 3 (c00). A seeded half (190 reads) gives 3 clusters (146 c01, 10 c00, 3 c01), same placements.
**Consensus placement** (identity x coverage; all 11 LRPAP1 loci on `pat`, the 9 mat sites on `mat`; `c01|c03` is one mat site because both lift to the same 20 kb):

| cluster (reads) | best pat locus | best mat site (nearest node) | next mat sites | O3 flag at 0.999 |
|---|---|---|---|---|
| 283, 1,492 bp | **c01 0.9987** | c00 0.9886 | p12 0.9880, c05 0.9873, c01\|c03 0.9839 | no (misses by 0.0003) |
| 4, 2,680 bp | **c01 0.9993** | c05 0.9888 | c04 0.9885, c02 0.9884, c01\|c03 0.9855 | **yes** |
| 21, 1,493 bp | c00 0.9980 | c00 0.9946 | | no (an allele of c00, 0.5% from the mother's) |
| 3, 1,496 bp | c00 0.9725 | c00 0.9732 | | no (matches neither) |

**Augmented reference** (reads aligned to the 4 consensus sequences; a read moves iff the consensus score is strictly higher than its maternal primary's). Half run, non-circular: 162 of the 191 held-out net reads move (84.8%), all 162 to the consensus of their own copy; median `de` 0.0115 -> 0.0027. The 2,494 reads already within the allele cutoff on the mother: 1,709 move (1,691 to the consensus of their own copy, 18 to another), median `de` 0.0040 -> 0.0007; 1,189 of 1,432 c00 reads and 505 of 530 c01 reads. Reads whose copy has no consensus (c02 to c05, p12/p14, p16, unplaced; 532 reads): 15 move (2.8%; c02 5, c05 10); p12/p14 (246), p16 (39), c03 (50), c04 (94) 0 moves.
**Bar.** A cluster flagged at the registered 0.999 with the mother below it: yes, the 4-read c01 cluster (0.9993 vs 0.9888). False moves <= 5% of the reads of the other copies: 2.8%, met. Verdict: **met, through the small cluster**; the 283-read cluster (96% of the c01 reads) has the same placement but 0.9987 and is not flagged at the registered bar.
**Not established / caveats.** (1) Truth is the paternal primary's locus (alignment-derived, not independent); p12 and p14 are tied and scored as one group. (2) That 1,691 placed reads move to a closer consensus is mostly the individual's own c00 and c01 sequences replacing the mother's: the move is toward the animal's own transcript, not a rescue of reads that were lost, so it is not counted as a gain. (3) CORRECTED after Amendment 9: the 283-read consensus aligns to pat c01 with identity 1.0000 (1,490 matches, 0 edits); the 0.9987 is coverage, 2 terminal consensus bases that do not align, not a consensus error. (4) The registered lift places c01 on the mat site it shares with c03 (222 mismatches in 20 kb); the consensus is 1.6% from that site (0.9839) and nearer to c00 (0.9886): the reads support c01 as the chr12 locus the mother's assembly lacks within the cutoff, which is evidence, not a test (single animal, one library).

---

### Amendment 9 (2026-10-08, before any polish number exists): an IsoCon-style correction step on the cluster consensus

**Why.** The consensus of a cluster reaches identity x coverage >= 0.999 against the erased copy in only 14 of 40 clusters on H and 65 of 212 on dev A (register drafts). In the LRPAP1 example the 283-read c01 consensus is at 0.9987, 2 bp short of the registered flag. IsoCon's second half is a correction step: align the cluster's reads to the candidate, test every non-candidate allele against the error rate, and change the candidate where the evidence is significant. The pool-first net so far stops at the abPOA consensus. This amendment tests whether that step closes the gap. It was thought of after seeing the LRPAP1 cluster (so LRPAP1 is not held out for it); the polish has not been run on A or H.
**Step (truth-free; reads and consensus only).** For each cluster: (1) align its reads to its consensus (`minimap2 -ax map-hifi --eqx`, primary records), count per consensus position the four bases and deletions, and per gap the inserted strings; (2) `e` = non-consensus observations / all observations of the cluster (a conservative estimate: it includes true variants); (3) an alternative allele with `k` of `n` reads is SIGNIFICANT iff `P(Binomial(n, e) >= k) < 0.05 / (3 L)` (Bonferroni over the three alternatives at each of the `L` positions; an insertion string is tested at its gap with the same correction); (4) a significant allele that is held by more reads than the consensus allele at that column REPLACES it (substitution, deletion of the base, insertion of the string); a significant allele held by fewer reads is only RECORDED (a minority variant, the IsoCon candidate for a split, not applied here); (5) repeat from (1) until a round changes nothing (at most 5 rounds). No other parameter.
**Scoring.** Unchanged (`fidelity.py`): the consensus (before and after) aligned with `minimap2 -c --cs -x splice:hq -uf -N 5` to the unmasked genome (`--cs` only adds a tag), best hit overlapping the erased interval of the cluster's majority family, identity x coverage >= 0.999. Dev A first, then held-out H; the code is committed after A and before H.
**Readings.** R1: clusters at >= 0.999 before and after (A, then H), clusters whose identity x coverage drops, total edit distance (NM) before and after. R2 (dev A only; the genome is used for this diagnostic, never by the polish): among the clusters still below 0.999 after the polish, the share of consensus-vs-genome edit columns that are significant minority-variant columns of the cluster's pileup (allele-like columns).
**Decision, before looking.** The correction step is supported iff on H the count at >= 0.999 rises from 14 to at least 19 of 40 (about a fifth of the 26 below the bar) AND at most 2 of the 40 clusters lose identity x coverage, AND the dev A count rises. If H rises by fewer than 5, the answer for the consensus-error explanation is NO; then, if R2 shows that at least half of the residual edit columns are significant minority-variant columns, the remaining shortfall is allelic (the reads carry two haplotypes, the assembly one) and allele-aware splitting is the next amendment (not run here); otherwise the shortfall is between the reads and the assembly and is not a consensus problem. The +5 bar is chosen here, with no polish data on A or H.
**Also reported (not judged):** the LRPAP1 consensus sequences after the polish, placed as in Amendment 8; run time.
**Not in scope:** splitting a cluster by its variants (needs a within-cluster truth); attribution (unchanged).

#### Amendment 9 result (2026-10-08, `bench/unmapped_rescue/{polish,run_polish}.py`; code committed after A and before H)

| | dev A (212 clusters) | held-out H (40 clusters) |
|---|---|---|
| identity x coverage >= 0.999, before -> after | 65 -> 70 (reproduces the stored 65) | 14 -> 15 (reproduces the stored 14) |
| clusters crossing the bar up / down | 8 / 3 | 1 / 0 |
| clusters whose identity x coverage drops / rises | 16 / 37 | **3** / 7 |
| total edit distance (NM) | 17,222 -> 15,942 (-7.4%) | 2,262 -> 2,094 (-7.4%) |
| corrections applied | 4,474 in 53 of 221 clusters (394 sub, 4,040 del, 40 ins) | 176 in 10 of 47 (16 sub, 159 del, 1 ins) |
| R2: residual edit columns that are significant minority-variant columns | 186 of 15,356 (1.2%) | 39 of 2,024 (1.9%) |

**Decision rule: NOT met.** H needed >= 19 of 40 (got 15) and <= 2 clusters losing identity x coverage (got 3); only the dev-A clause (count rises) holds. The consensus-error explanation of the shortfall is **refuted**. R2 is far below half, so the residual is not allelic either; per the rule the shortfall lies between the reads and the assembly, not in the consensus.
**What the corrections were.** 90% are deletions, concentrated in a few clusters whose consensus carried an exon most reads skip (e.g. 4,912 -> 2,882 bp, 2,030 deleted positions): a mixed-isoform cluster trimmed to its plurality isoform, not base-level repair. Base-level substitutions were 394 of 4,474 on A and 16 of 176 on H.
**Exploratory, on dev A only (not a result).** The metric is dominated by something polishing cannot touch: identity alone reaches 0.999 in 95 of 212 clusters but identity x coverage in 65, so 30 clusters are lost to coverage (unaligned consensus ends) alone; the median identity x coverage is 0.9976; 22 clusters have more than 200 edits against the assembled copy (the cluster is not a single copy end to end); NM is 0 in more than a tenth of the clusters. The registered 0.999 bar measures end-to-end agreement with the annotated copy more than consensus accuracy.
**LRPAP1 (reported, not judged).** The polish applied 0 corrections to the 4 consensus sequences; placements are identical to Amendment 8. Checked against the paternal genome: the 283-read consensus aligns to pat c01 over 1,490 of its 1,492 bases with 0 edits (identity 1.0000), so the 0.9987 is 2 unaligned terminal bases. Its 21 "minority variants" are single-base deletions held by 8 to 11 of about 275 reads each (HiFi dropout), not alleles.
**Stage 2 (allele-aware splitting) is not warranted by this amendment's rule** (R2 1.2% and 1.9%).

---

### Amendment 10 (2026-10-08, before any synthetic number exists): a supervised synthetic world, erase one copy, recover it

**Why.** Two things the real beds cannot answer: (i) what the unaligned consensus ends (1 to 4 bp in 21 of the 30 coverage-limited dev-A clusters, 5' extensions of 18 to 1,040 bp in 4; no polyA) are, because the true transcript is unknown there; (ii) how recovery depends on how close the erased copy is to a surviving copy, because the real families mix every divergence. A synthetic world gives exact truth: every read's copy, every copy's transcript and genomic sequence.
**World** (`bench/unmapped_rescue/synth_world.py`, built with `bench/famsim`; seed 20261008; nothing tuned afterwards). 30 families, 6 for each divergence class D in {0.5%, 1%, 2%, 4%, 8%}. Per family a random synthetic gene (5 to 10 exons, internal exons 90 to 350 bp, first 100 to 250, last 250 to 800, introns 300 to 2,500 bp, GT...AG), three copies on independent random background: A (the template, present), B (present, SNP rate 1.5 D), E (ERASED: the contig is omitted from the reference, SNP rate D; its nearest survivor is A at about D). 80 reads per copy (HiFi substitutions 0.001, short indels 0.0003, end jitter 30 bp at both ends; famsim's read model, which has no polyA, adapter or 5' degradation). Read-end variants on the same world: E0 as above; E1 = E0 + a 3' polyA tail of 20 to 30 A on every read; E2 = E0 + 5' truncation (30% of reads lose up to 30% of their length).
**Pipeline (the frozen one, truth-free):** reads aligned to the reference with the study's command; the net = unmapped or divergence > 0.00958; seed-round clustering with the frozen edge rule, components >= 3; abPOA consensus; dc-megablast cover-score attribution against the surviving copies (margin 1.10); realignment of all reads to the consensus sequences (Amendment 7 rule).
**Readings and predictions (bars set now):**
- S1 (the simulation behaves): for erased copies with D >= 2%, >= 95% of their reads are in the net; for D = 0.5%, < 20% (designed-undetectable: the cutoff is about 1%).
- S2 (consensus): for D >= 2%, a cluster of >= 3 reads exists for >= 90% of the erased copies, purity 1.000 by true copy, and its consensus aligns to the TRUE transcript of the erased copy with identity >= 0.999 over the aligned part in >= 90% of those clusters.
- S3 (the registered metric): identity x coverage >= 0.999 of the consensus against the true genomic copy (the metric used on the real beds), same clusters. Reported next to S2. If S3 << S2 the ends are the cause.
- S4 (ends, descriptive): per end, the consensus overhang beyond the true transcript and the true bases it lacks, and the aligner's soft clip, per read-end variant. Prediction: with jitter only, the consensus ends fall inside the read-end envelope (within 30 bp of the true ends); polyA (E1) appears as a 3' overhang of about 25 bp; the 1 to 4 bp clips of the real beds are NOT reproduced by E0.
- S5 (attribution): net reads of D >= 2% erased copies attributed to the correct family >= 90%, wrong <= 1%.
- S6 (augmented reference): erased-copy reads in the net land on their own consensus >= 90%; reads of surviving copies move <= 1%.
- S7 (the divergence curve, descriptive): S1, S2, S5 and S6 by D.
**Decision.** The supervised check "the pipeline recovers an erased copy when it is recoverable" passes iff S1, S2, S5 and S6 hold on E0. The ends question is answered by S3 vs S2 and S4: if E0 reproduces the real-bed signature (S3 << S2, clips of 1 to 4 bp) the mechanism is in the pipeline (aligner or consensus end handling); if only E1 or E2 reproduces it, it is in the reads; if none does, the real-bed ends come from something the simulation lacks (the individual's variation against the assembly, UTR structure), and the next step is read-end inspection on real data.

#### Amendment 10 result (synthetic world, 2026-10-08; `synth_world.py`, `run_synth.py`, `ends.py`; seed 20261008; famsim `verify` 28/28 PASS on the families checked)

30 families x 3 copies (A, B present; E erased), 80 reads per copy, 7,200 reads; E0 = jitter only. The net is every read unmapped or with divergence > 0.00958.

| D (E vs nearest survivor) | E reads in the net (E0) | E copies with a pure cluster | consensus vs TRUE transcript, core identity >= 0.999 | identity x coverage >= 0.999 vs true genome copy (E0 / E1 polyA / E2 trunc) | attributed to the right family, wrong | E reads landing on own consensus | survivor reads moved, AS rule / divergence rule |
|---|---|---|---|---|---|---|---|
| 0.5% | 12.1% (E2 14.8%) | 2 of 6 | 2 of 2 | 2/2, 0/2, 2/2 | 57 of 58, 0 | 57 of 58 | 306 of 960 / 0 |
| 1% | 99.2% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 476 of 476, 0 | 476 | 480 of 960 / 0 |
| 2% | 100% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 480 of 480, 0 | 480 | 289 of 960 / 0 |
| 4% | 100% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 480 of 480, 0 | 480 | 0 / 0 |
| 8% | 100% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 480 of 480, 0 | 480 | 0 / 0 |

**S1 holds** (D >= 2%: 100% in the net; D = 0.5%: 12.1%, designed-undetectable; at 1% the cutoff is crossed, 99.2%). **S2 holds** (every D >= 1% copy has a pure cluster and its consensus is >= 0.999 identical to the true transcript over the aligned part; the 0.5% class has clusters for only 2 copies because only 58 of its 480 reads are in the net). **S5 holds** (0 wrong). **S6 FAILS under the registered score rule**: reads of the SURVIVING copies move onto the erased copy's consensus in 30 to 50% of cases for D <= 2% (bar <= 1%); they are all A reads (de 0.001 to 0.003 on the genome, AS 1221 to 1248) that score higher on the unspliced consensus (AS 1279, de 0.0138): a spliced genome alignment pays about 8 points per junction that an unspliced transcript alignment does not, so raw scores are not comparable. **Registered decision (S1, S2, S5, S6 all hold on E0): NOT met, because of S6 only;** recovery itself (S1, S2, S5) holds. The failure is in the move rule of Amendment 7; real bed H (D mostly >= 2%) did not show it (0 of 41,727), but a family whose erased copy is within about 2% of a survivor would.
**S3/S4, the ends question:** E0 does not reproduce the real-bed signature: identity x coverage equals the core identity (6/6), the aligner clips 0 bases (median), and the consensus is SHORTER than the true transcript (5' lacks 0 to 23 bases, 3' lacks 0 to 7: the read-end jitter), never longer. E2 (5' truncation) is the same. E1 (polyA) fails the metric completely (0 of 30 families; the 3' tail is 13 to 30 bases of the consensus, clipped by the aligner at a median of 24 to 28) which the real beds do not show (no A-rich clips). **None of E0, E1, E2 reproduces the real beds, so by the registered rule the real-bed ends come from something the simulation lacks; real-data inspection is the next step** (Amendment 11).

---

### Amendment 11 (2026-10-08, AFTER the Amendment 10 numbers; exploratory, post hoc, labelled as such)

Two things the synthetic world exposed, followed up on the real beds. Neither is a registered test; both are descriptive, dev A first. H was already used by the earlier amendments, so H numbers here are not held-out.
**(a) A divergence move rule.** A read moves to a consensus iff the consensus alignment's `de` is strictly lower than its genome primary's (unmapped reads: as before). Synthetic world: survivors moved 0 of 960 in every class and variant (AS rule: up to 480), erased-copy reads in the net still all land on their own consensus, and in the 0.5% class (copy hidden inside the allele cutoff) 103 of 422 absorbed erased-copy reads move to their own consensus where one exists. Real bed H, half run: held-out unmapped reads 2,505 of 2,657 (identical to the score rule), absorbed 2,122 of 11,973 (score rule 2,250), survivors 0 of 41,727: H's result does not depend on the rule.
**(b) The unaligned consensus ends of the real beds are an untemplated 5' run of G.** In consensus (transcript) orientation, on dev A 129 of 212 consensus sequences have a 1 to 20 bp unaligned 5' end, 70% G, 102 of them pure G; lengths 1 bp (56), 2 bp (44), 3 or more (29). None of 127 clips of 4 bp or less equals the adjacent genomic sequence (a random same-size window matches 24 of 127), so the aligner is not clipping genomic bases: the consensus carries bases that are not in the genome. 3' clips are rare (20) and not composition-biased. This is the template-switching / cap artifact (a few untemplated G at the 5' end of the cDNA); the earlier read-ends note (`project_read_proven_ends`) recorded no such signal for gorilla, which this contradicts for these consensus sequences and needs reconciling. Excluding a pure 5' G run of at most 3 bp from the coverage denominator turns identity x coverage >= 0.999 from 65 to 80 of 212 clusters (dev A) and from 14 to 17 of 40 (H): a metric that ignores the artifact would score the consensus higher, not a change to any registered result. The remaining shortfall is internal identity (identity alone >= 0.999 in 95 of 212) and 5' extensions of 18 to 1,040 bp (4 of 30 coverage-limited clusters).
**Next, proposed (not run):** add the untemplated 5' G run to the read simulator so the supervised world reproduces the real signature, then test trimming it from the consensus with exact truth.

---

### Amendment 12 (2026-10-08, before any number from the fresh world): adopt the divergence move rule in the augmented reference

**Decision.** The augmented-reference step (Amendments 7, 8, 10) moves a read to a consensus iff the consensus alignment's divergence `de` is strictly lower than its genome primary's (an unmapped read: consensus at identity >= 0.98 and query cover >= 0.80, as before). It replaces the raw-score rule, which compares a spliced genome alignment with an unspliced transcript alignment (about 8 points per junction in favour of the transcript) and moved 30 to 50% of surviving-copy reads in the synthetic world when the erased copy was within 2% of a survivor. The score rule stays available as an explicit option (`rule="score"`) so earlier results can be reproduced. The rule was found on the synthetic world of seed 20261008 (Amendment 10), so that world is no longer held out for it.
**Test (fresh world, never looked at).** `synth_world.py build --seed 20261009` with eight divergence classes, four families each (32 families): D in {0.5%, 0.75%, 1%, 1.5%, 2%, 3%, 4%, 8%}, three of which (0.75%, 1.5%, 3%) were not in the first world. Everything else as in Amendment 10 (E0 read model, frozen pipeline).
**Bars (set now).** (1) Survivor reads moved onto an erased-copy consensus <= 1% in every class (the score rule gave up to 50%). (2) Erased-copy reads in the net land on their own family's consensus >= 90% in every class with a pure cluster for at least 3 copies. (3) Recovery is unchanged by the rule: copies with a pure cluster and attribution results identical to the score-rule run on the same world (the rule only touches the move step). Reported, not judged: erased-copy reads absorbed inside the allele cutoff (not in the net) that move to their own consensus, by class.
**Also re-run under the new default (reported):** bed H half run (previous: held-out unmapped 2,505 of 2,657, absorbed 2,122, survivors 0 with the divergence rule) and the LRPAP1 example (Amendment 8).
**Decision rule.** The adoption stands iff bars (1) to (3) hold. If (1) fails in any class, the rule goes back to explicit-option status and the failure is reported.

#### Amendment 12 result (fresh world seed 20261009, 32 families; 2026-10-08)

Net 2,015 of 2,560 erased-copy reads, 0 survivor reads; 26 clusters, 26 attributed.

| D | E reads in net | E copies with a pure cluster | survivor reads moved, divergence rule (score rule) | E net reads on own consensus |
|---|---|---|---|---|
| 0.5% | 0% | 0 of 4 | 0 of 640 (0) | |
| 0.75% | 54.7% | 3 of 4 | **0** (320) | 175 of 175 |
| 1% | 75.0% | 3 of 4 | **0** (240) | 240 of 240 |
| 1.5% | 100% | 4 of 4 | **0** (241) | 320 of 320 |
| 2% | 100% | 4 of 4 | **0** (71) | 320 of 320 |
| 3%, 4%, 8% | 100% | 4 of 4 | 0 (0) | 320 of 320 |

**Bars: (1) survivors moved 0 of 640 in every class (<= 1%): met. (2) erased-copy net reads on own consensus 100% in every class with >= 3 pure clusters: met. (3) recovery unchanged: the move rule acts after clustering and attribution, which are identical to the score-rule run: met.** The adoption stands; `augment.new_primary` and `move_metrics` now default to the divergence rule (`rule="score"` reproduces Amendments 7, 8, 10).
The detection edge sits where the cutoff says: 0% of the erased copy's reads leave it at 0.5%, 55% at 0.75%, 75% at 1%, all from 1.5%. Reported: erased-copy reads inside the cutoff that move to their own consensus: 65 of 145 at 0.75% (where a consensus exists).
**Re-runs under the new default.** Bed H half run: held-out unmapped 2,505 of 2,657 (94.3%), absorbed 2,122 of 11,973, survivors 0 of 41,727, background onto a family cluster 0 of 959: the H conclusions do not depend on the rule. LRPAP1 (Amendment 8): held-out net reads 172 of 191 move (90.1%; 162 under the score rule), already-placed reads 1,121 move (1,709 under the score rule; 1,118 to the right copy), other-copy reads 3 of 532 (0.6%; 15 under the score rule).

---

### Amendment 13 (2026-10-08, before any number): the untemplated 5' G run, simulated and trimmed

**What is known (Amendment 11, dev A).** 129 of 212 consensus sequences have an unaligned 5' end; 102 are pure G, 1 to 2 bp, none matches the adjacent genome. From the reads of the 53 dev-A clusters (>= 10 reads) whose consensus carries a 1 to 3 bp 5' G run: the median cluster has 88% of its reads starting with that run (quartiles 74%, 94%); over all their reads the leading G run is 0 bp in 5%, 1 bp 21%, 2 bp 41%, 3 bp 23%, 4 bp 7%, 5 or more 3%. Reads are in transcript orientation (the consensus is, and its 3' end is not polyA).
**Simulation (read variants on a fresh world, seed 20261010, same 30-family design as Amendment 10).** E4 = reads with end jitter 3 bp (a dominant start site, unlike E0's 30; no other change). E3 = E4 plus an untemplated 5' G run on 95% of reads, run length drawn from the measured distribution above (1: 0.22, 2: 0.43, 3: 0.24, 4: 0.07, 5: 0.03, 6: 0.01 of the reads that carry it).
**Trimming rule T (truth-free; uses what the pipeline already has).** For a cluster attributed to a family: take the surviving copy of that family with the largest union of dc-megablast HSP query spans (the cover score's best copy); the 5' PREFIX of the consensus is its bases before the first HSP's query start on that copy. If the prefix is 1 to 3 bases and all G, remove it from the consensus; otherwise change nothing. Unattributed clusters are not touched. No other parameter.
**Readings and bars (set now), classes D >= 1% (30-D classes with a pure cluster):**
- (1) E3 reproduces the real signature: identity x coverage >= 0.999 against the TRUE genome copy in <= 50% of E3 clusters (real beds: 31% dev A, 35% H) while >= 90% in E4 (the control). If not, the simulation does not capture the artifact and the rest is not interpreted.
- (2) T repairs it: after T, identity x coverage >= 0.999 in >= 90% of E3 clusters, and the 5' offset of the consensus against the TRUE transcript (ends.py) within [-4, +1] in >= 90%.
- (3) T does no harm: in E4 (no artifact) T trims at most 1 cluster; in E3 T trims a base that belongs to the true transcript (5' offset after T below -4) in at most 1 cluster.
- (4) Reported: T's decisions on E3 by D class (trimmed, not trimmed because abstained, not trimmed because prefix too long or not G).
**Real data (descriptive, dev A first; H is not held out for diagnostics).** Apply T with the real attribution and the real HSPs: how many consensus sequences it trims, agreement with the genome-based finding (a pure-G 1 to 3 bp 5' clip against the erased copy): precision, recall; and identity x coverage >= 0.999 after T by excluding the trimmed prefix from the denominator (no re-alignment needed when the prefix is within the unaligned clip).
**Decision.** T is adopted as a consensus clean-up step iff (1) holds and (2), (3) hold. If (1) fails, no conclusion about T; if (2) or (3) fails, T is reported as not adopted with the failure.

#### Amendment 13 result (fresh world seed 20261010; 24 erased copies with a pure cluster at D >= 1%; 2026-10-08)

| | E4 (jitter 3, no artifact) | E3 (E4 + 5' G run) untrimmed | E3 after T |
|---|---|---|---|
| identity x coverage >= 0.999 vs true genome copy | 24 of 24 | **14 of 24 (58%)** | **24 of 24** |
| 5' offset vs true transcript within [-4, +1] | 24 of 24 | 21 of 24 | 24 of 24 |
| consensus with an unaligned 5' clip | 0 of 24 | 20 of 24 (lengths 2: 12, 1: 7, 3: 1) | |
| clusters trimmed by T / over-trimmed (offset below -4) | 0 / 0 | | 21 / 0 |

**Bar (1) NOT met (58% > 50%):** E3 produces the 5' clip itself (20 of 24 consensus sequences, real dev A 129 of 212) but, with an error-free consensus, a 1 to 2 bp clip on a 2.5 kb transcript scores 0.9996 or 0.9992 and still passes 0.999; only runs of 3 or a short transcript fail. The real beds fall to 31 to 35% because the clip adds to other small deviations. Bars (2) and (3) hold. **By the registered decision (needs (1)), T is not adopted.**
**On the real beds T does almost nothing (descriptive):** dev A: trimmed 4 of 212 (3 agree with the genome-based pure-G clip), prefix longer than 3 in 115, abstained 87; H: trimmed 0 of 40 (prefix too long 25, abstained 13). The dc-megablast HSP does not start at the first consensus base on a real survivor (a short or divergent first exon starts the HSP tens of bases in), so the HSP prefix is the wrong measurement outside the synthetic world, where the first exons are identical. identity x coverage >= 0.999 unchanged (65 and 14).

---

### Amendment 14 (2026-10-08, after the Amendment 13 numbers, before any T2 number): rule T2, aligner-based prefix

**Change.** The 5' prefix is measured by aligning the consensus (minimap2 `-c -x splice:hq -uf -N 5`, the command used for the fidelity metric) to the genomic spans of the surviving copies of its attributed family and taking the consensus bases before the alignment start (PAF query start) of the best record by matches; the decision step is unchanged (1 to 3 bases, all G, else no change; abstained clusters untouched). Nothing else changes.
**Test.** (a) Fresh synthetic world, seed 20261011 (same design, E3 and E4): bars (2) and (3) of Amendment 13, and bar (1') replacing (1): E3 shows an unaligned 5' clip in >= 50% of the consensus sequences and E4 in 0 (the property the real beds have: 61% dev A). (b) Real dev A, descriptive but with a bar for the decision: among attributed clusters that have a genome-based pure-G 1 to 3 bp 5' clip (against the erased copy), T2 recovers >= 80% of them, and >= 95% of T2's trims coincide with such a clip. H and the identity x coverage count after T2 are reported, not judged.
**Decision.** T2 is adopted as a consensus clean-up iff (1'), (2), (3) hold on the synthetic world AND the dev-A recall and precision bars hold. Otherwise it is reported as not adopted with the failure.

#### Amendment 14 result (2026-10-08; synthetic world seed 20261011, real dev A and H)

**Synthetic (23 erased copies with a pure cluster at D >= 1%).** E4 (no artifact): identity x coverage >= 0.999 in 23 of 23, T2 trims 0, no over-trim (bar (3)). E3: an unaligned 5' clip in 22 of 25 consensus sequences (E4: 0 of 25), bar (1') met. Untrimmed E3: 12 of 23 (52%) reach 0.999; **after T2: 21 of 23 (91%)**, 5' offset within [-4, +1] in 22 of 23 (96%), over-trim 0, 18 clusters trimmed (the 2 failures: one prefix longer than 3, one not pure G). Bars (1'), (2), (3) met.
**Real dev A (the decision bar): not met.** T2 acts on 3 of 212 consensus sequences: `no_hsp` 122 (no alignment to any surviving copy of the attributed family), `abstained` 87, `too_long` 3. Of the 60 attributed clusters that have a genome-based pure-G 5' clip, T2 trims 0 (recall 0%, bar 80%). H: 0 of 40 (13 attributed clusters with a clip, 0 trimmed). The cause is alignment sensitivity: `minimap2 -x splice:hq` against the survivors yields any record for only 44 of 221 consensus sequences (median block 280 bp), because real survivors are 3 to 15% diverged paralogs, not the near-identical copies of the synthetic world. identity x coverage >= 0.999 unchanged (65, 14).
**Decision: T2 is NOT adopted** (the synthetic bars hold, the real-data recall bar fails). What the synthetic world establishes: when the survivors align at the first exon, the pure-G prefix is identifiable and trimming it restores the registered metric with no harm. What it does not: a way to identify the artifact on real families whose survivors are too divergent to anchor the 5' end. The artifact is a property of the library (a 5' soft-clipped pure-G run on reads that align cleanly elsewhere), which is not used by T or T2.

---

### Amendment 15 (2026-10-08, before any T3 number): rule T3, the artifact read from the library

**What the library shows (design check on the 41,727 surviving-copy reads of the real library, `linktest/R.bam`; looked at before this amendment).** Of 39,351 cleanly aligned survivor reads (primary, MAPQ >= 10, divergence <= 0.0096), 28,581 (72.6%) have a 5' soft clip that is a pure G run of 1 to 3 bases (1 bp 12,770; 2 bp 14,843; 3 bp 968); the same clip made of A, T or C: 62, 12, 5. The artifact is a property of the library, readable from reads that align perfectly elsewhere.
**Rule T3.** (1) LIBRARY GATE: from the cleanly aligned reads of the surviving copies count the 5' soft clips (read orientation) that are a pure single base of length 1 to 3; the artifact is declared present iff the G count has P(Binomial(n, 0.5) >= G count) < 1e-6, n = the pure-A/C/G/T clip count (needs n >= 20), i.e. G is not just one of four equally likely bases. (2) If present, remove the maximal leading G run of every consensus sequence if its length is 1 to 3 (a longer run is left alone); no reference, no attribution. If absent, nothing is trimmed.
**Test.** (a) Fresh synthetic world, seed 20261012 (E3 and E4; the gate is computed from the survivor reads of that world's own reference alignment). Bars: gate TRUE on E3 and FALSE on E4; E4: 0 consensus sequences trimmed; E3: identity x coverage >= 0.999 vs the true genome copy in >= 90% of the D >= 1% erased-copy clusters and the 5' offset against the true transcript within [-4, +1] in >= 90%. (b) Real dev A (decision bars, genome-based finding as reference): gate TRUE; of the consensus sequences with a pure-G 1 to 3 bp 5' clip against the erased copy, T3 trims >= 95%; among T3's trims, the trimmed run exceeds the unaligned clip by at most 1 base in >= 90% (templated G's are not lost); identity x coverage >= 0.999 after T3 (excluding the trimmed bases from the denominator, trimmed bases that the aligner had covered subtracted from the aligned length) >= 75 of 212. H reported, not judged (descriptive prediction 16 or more of 40).
**Decision.** T3 is adopted as the consensus clean-up iff the gate behaves (TRUE on E3 and real, FALSE on E4) and every bar in (a) and (b) holds; otherwise not adopted, with the failure.

#### Amendment 15 result (2026-10-08; synthetic world seed 20261012; real library; dev A and H)

**Gate.** Real library: OPEN (G 28,614 of 28,694 pure 5' clips; A 63, T 12, C 5; 39,403 clean survivor reads; p = 0). Synthetic E3: OPEN (G 3,829, others 0); E4: closed (A 3, C 5, G 2, T 1; p = 1). The gate behaves.
**Synthetic (D >= 1%, 24 erased copies with a pure cluster).** E4: 0 trimmed, 24 of 24 at 0.999 untouched. E3: identity x coverage >= 0.999 in **12 of 24 -> 22 of 24 (92%)**, 5' offset within [-4, +1] in 14 -> **23 of 24 (96%)**, over-trim 0; the two untrimmed have a leading run longer than 3. Bars (a) met.
**Real dev A (decision bars): two of three FAIL.** T3 trims 197 of 212 consensus sequences (too_long 13, no run 2). Recall among the 100 consensus sequences with a genome-based pure-G 1 to 3 bp clip: **92 of 100 (92%; bar 95%) FAILS** (the other 8 have a leading run of more than 3). Templated-sequence check: the trimmed run exceeds the unaligned genome clip by at most 1 base in **132 of 197 (67%; bar 90%) FAILS** (excess 0: 39, 1: 93, 2: 51, 3: 14). identity x coverage >= 0.999: 65 -> **82 of 212** (bar >= 75) met. H (descriptive): 14 -> 18 of 40; trims 37 of 40, 19 of the 20 with a genome clip.
**Decision: T3 is NOT adopted** (two of the three registered real-data bars fail). What this shows: the artifact is a library property that a gate can read perfectly and that a leading-G trim removes, restoring the registered metric on synthetic worlds (and on the real beds, 65 -> 82 and 14 -> 18). What it cannot do is separate artifact G from templated G: real 5' ends are GC-rich, and 105 of the 197 trimmed consensus sequences have no unaligned G clip at all, so their leading G run was aligned to the erased copy (templated, or artifact G that matched the G-rich upstream flank by chance, which the alignment cannot tell apart). The registered over-trim bar uses the unaligned clip as the artifact's length, which therefore under-counts it; on real data the bar was mis-specified for what can be observed, but it was registered and it failed.

---

### Amendment 16 (2026-10-08, before any T4 number): rule T4, per-cluster correction of the 5' G run from the reads' own leading-G lengths

**Why.** T3 trims every leading G run of 1 to 3 and cannot tell the artifact's G's from the transcript's own (real 5' ends are GC-rich). The reads can: every read of a cluster carries the same templated run of j G's, then an artifact run a_i whose length distribution is a property of the library. So a read's leading G run is r_i = j + a_i (read starts are within a base or two of each other; the model is deliberately simple).
**Rule T4 (gate as in Amendment 15, nothing else added).** (1) Library artifact length distribution a(l), l = 0..8, from the cleanly aligned survivor reads: a(l) = (reads whose 5' clip is a pure G run of l bases) / (clean reads), a(0) = reads with no clip / clean reads; a floor of 1e-4 for lengths never seen. (2) For each cluster with >= 5 reads: r_i = the leading G run of each read's own sequence (no alignment, no CIGAR); j_hat = argmax over j in 0..6 of the sum of log a(r_i - j), a term with r_i < j taking 1e-6; ties go to the larger j. (3) The consensus has a leading G run of k; remove the first max(0, k - j_hat) bases (all G). Gate closed, or a cluster of fewer than 5 reads: no change.
**Test world (fresh, seed 20261013; 30 families, same classes).** Templates whose first exon begins with j = 0, 1, 2 or 3 G's (forced, equally many families) to reproduce the GC-rich 5' ends; E5 = end jitter 1 + the measured artifact (95% of reads, mode 2); E6 = the same without the artifact (control). Truth: s* = the modal read start (chain offset) among the cluster's reads; the correct consensus is the true transcript from s* on, so the 5' offset against the true transcript should be -s*; err = (offset + s*), 0 exact, > 0 artifact left in, < 0 templated bases removed. T3 (trim all) and the untrimmed consensus are scored on the same clusters.
**Bars (set now), D >= 1% erased-copy clusters:** (1) gate open on E5, closed on E6 (E6: nothing trimmed). (2) T4 err = 0 in >= 80% and |err| <= 1 in >= 95%. (3) T4 removes templated bases (err < 0) in <= 10% and leaves artifact (err > 0) in <= 10%. (4) T4 has strictly fewer err < 0 clusters than T3. (5) identity x coverage >= 0.999 against the true genome copy in >= 90% of clusters after T4.
**Real dev A (one bar, rest descriptive):** among consensus sequences with a genome-based pure-G 5' clip of g bases, T4 trims at least g in >= 90% (it does not leave unaligned artifact G's); identity x coverage >= 0.999 after T4 >= 75 of 212; reported: distribution of j_hat and of the trimmed amount against T3's. H descriptive.
**Decision.** T4 is adopted as the consensus clean-up iff bars (1) to (5) and the dev-A bar hold; otherwise not adopted, with the failure.

#### Amendment 16 result (2026-10-08; synthetic world seed 20261013 with forced leading G; real dev A and H)

**Scoring check first.** The control E6 (no artifact, same templates) is exact in 24 of 24 clusters under all three treatments, so the err measure (5' offset against the true transcript + the modal true read start) is sound. (A first scoring run used the E0 read-start table by mistake, gave nonsense errors up to 30, and was discarded before any reading; the per-variant truth tables fixed it, the world is byte-identical.)
**Synthetic E5 (artifact + 0 to 3 templated leading G), D >= 1%, 24 erased-copy clusters:**

| | untrimmed | T3 (trim all, 1 to 3) | **T4** |
|---|---|---|---|
| err = 0 (exact) | 0 | 7 | **14 (58%)** |
| \|err\| <= 1 | 4 | 16 | **24 (100%)** |
| templated bases removed (err < 0) | 0 | 7 | **1 (4%)** |
| artifact G left (err > 0) | 24 | 10 | **9 (38%)** |
| identity x coverage >= 0.999 | 16 | 20 | **24 (100%)** |

j_hat equals the forced templated count or one less in every cluster (6, 6, 6, 7 clusters for j = 0 to 3). The 9 clusters with a G left all have the modal read start at transcript position 1 (jitter 1: half the reads start at 0, half at 1), where the templated count at the consensus start is one lower than the global estimate: a +1 ambiguity of the start itself.
**Bars: (1) met** (gate open on E5, closed on E6, nothing trimmed there). **(2) FAILS:** exact 58% (bar 80%); |err| <= 1 in 100% (bar 95%) met. **(3) FAILS:** templated removed 4% (met), artifact left 38% (bar 10%). **(4) met** (1 against T3's 7). **(5) met** (100%).
**Real dev A (bar: T4 removes at least the unaligned genome clip in >= 90%): FAILS, 64 of 100 (36% under-trimmed).** Of the 100 consensus sequences with a genome G clip, 77 come from clusters of >= 5 reads (T4 applicable): 64 ok, 13 short by one base; the other 23 are clusters of fewer than 5 reads, where T4 does nothing by rule. identity x coverage >= 0.999: 65 -> **80** of 212 (bar >= 75, met). j_hat: 0 in 98 clusters, 1 in 51, 2 in 4 (no estimate for 59 small clusters). T3 would remove more than T4 in 104 clusters (1 base: 65, 2: 30, 3: 9). H (descriptive): 14 -> 17 of 40; 6 of 20 clipped consensus sequences short (4 of them clusters < 5 reads).
**Decision: T4 is NOT adopted by the registered rule** (bars (2), (3) and the dev-A bar fail). In its favour and descriptively: it is the only treatment with |err| <= 1 everywhere, it removes templated bases 6 times less often than T3, it restores the metric in all synthetic clusters, and on dev A its shortfall is the 1-base start ambiguity (13 of 77) plus clusters too small for the estimate (23 of 100). The bars of (2) and (3) asked for exactness that the 1-base jitter of the read starts rules out.

---

### Amendment 17 (2026-10-08, before any gate-aware flag number): the registered identity x coverage ignores a leading G run when the library gate is open

**Change (score side only; the consensus is never touched).** The registered metric (identity x coverage >= 0.999 of the consensus against a genome copy, which is also the O3 flag's recovery hit) is ident x (aligned query bases / consensus length). Gate-aware version: if the library gate of Amendment 15 is open AND the consensus' unaligned 5' end (the bases before the alignment start, consensus = transcript orientation) is 1 to 3 bases, all G, those bases are removed from the denominator; in every other case the metric is unchanged. A 5' end that is not pure G, longer than 3, or any gate-closed library is scored exactly as before. 3' ends are never touched.
**Tests.** (a) Synthetic, fresh data not used for this rule: the worlds of seeds 20261012 (E3 artifact, E4 control) and 20261013 (E5 artifact + 0 to 3 templated leading G, E6 control), D >= 1% erased-copy clusters. Bars: gate-aware >= 0.999 in >= 90% of the clusters of E3 and of E5 (registered metric: 12 of 24 and 16 of 24); control worlds (gate closed) identical to the registered metric in every cluster; **no new false passes**: every cluster that passes only under the gate-aware metric has core identity >= 0.999 against the TRUE transcript (0 exceptions). (b) Unit test: a non-G or longer-than-3 5' end is scored as before. (c) Real dev A and H, reported not judged (the 80 and 17 of Amendment 11 were computed with the same exclusion from the genome clip, so these are not blind): gate-aware count. (d) LRPAP1 (Amendment 8): the O3 flag of each of the four consensus sequences, registered and gate-aware, from the existing alignments to `pat` and `mat` (a cluster's flag = best pat locus >= 0.999 and best mat site < 0.999).
**Decision.** Adopted as the flag/fidelity metric iff every bar in (a) and (b) holds.

#### Amendment 17 result (2026-10-08)

**(a) Synthetic, D >= 1% erased-copy clusters (24 per variant):** E3 (artifact): registered 12 -> **gate-aware 23 (96%)**; E5 (artifact + 0 to 3 templated leading G): 16 -> **24 (100%)**; controls with the gate closed (E4, E6): identical in 24 of 24 clusters (24 of 24 at 0.999 both ways); **new passes without core identity >= 0.999 against the true transcript: 0** in all four. **(b)** unit tests: a non-G or longer-than-3 5' end, a closed gate, and a consensus with a real mismatch (99.5% identity) are scored as before. Bars (a), (b) met: **the gate-aware metric is adopted** (`flagmetric.gate_aware`; Python-side scoring and reporting only; the Rust `o3_candidates` flag is not changed).
**(c) Real beds (not blind: the same exclusion gave 80 and 17 in Amendment 11):** dev A 65 -> 80 of 212, H 14 -> 17 of 40.
**(d) LRPAP1 O3 flag (best `pat` locus >= 0.999 and best `mat` site < 0.999), registered -> gate-aware:** 283-read c01 consensus: pat c01 0.9987 -> **1.0000**, mat best c00 0.9886 -> 0.9899: **not flagged -> flagged**. 4-read c01 cluster: 0.9993 -> 1.0000 vs mat 0.9895, flagged both ways. 21-read c00 cluster: pat c00 0.9980 -> 0.9993 vs mat c00 0.9946 -> 0.9959: **not flagged -> flagged**; this is the c00 allele of Amendment 8, 0.4% from the mother's, so the flag as defined does not separate an allele from an absent copy (it never did; the registered bar is identity, not the 0.958% allele cutoff). 3-read cluster: its 5' end (GAG) is not pure G, unchanged, not flagged.

---

### Amendment 18 (2026-10-08, before any partition number): IsoCon-style partition of a cluster into competing candidate transcripts

**Why.** The frozen pipeline reports ONE consensus per cluster. A cluster that mixes two isoforms (exon skipping) or two near-identical copies (< 0.958% apart, joined by the edge rule) is reported as the majority isoform / a column-wise blend, and the other member is lost. Amendment 9 tested base-level correction (refuted) and the allele-like columns (1 to 2%), not the separation of isoform- or copy-scale members. IsoCon's remaining idea is exactly that: candidates compete for reads and a candidate stands only if statistically supported.
**Algorithm PART (no reference, no attribution).** For a cluster of >= 6 reads (at most 400, seeded sample for discovery): (1) its consensus C (abPOA of the 100 longest reads, as frozen); (2) align the reads to C (`minimap2 -ax map-hifi --eqx`), take the significant variant columns V of the Amendment 9 test (substitution, deletion, insertion alleles at Bonferroni 0.05; both the alleles the consensus carries as minority and those it lacks), dropping the first and last 15 consensus columns (read ends carry the jitter and the 5' artifact); (3) for every pair of variants in different columns, a one-sided hypergeometric (Fisher) test of positive association over the reads covering both, edge iff p < 0.05 / (number of pairs) and both carried by >= 3 reads; (4) a BLOCK = a connected component of >= 2 variants (a lone significant column, like a sequencing hotspot, is not a block); a read carries a block if it carries at least half of the block's variants that it covers; (5) split on the block whose smaller side is largest, if both sides have >= 3 reads; each side gets its own abPOA consensus and the procedure repeats (depth <= 6) until no block remains. The leaves, with their reads, are the candidates. Parameters: alpha 0.05, 15-column end window, 3 reads (the frozen minimum cluster), 6 reads to split, depth 6.
**Synthetic world (seed 20261014; `synth_world.py partition_specs`; 12 families per scenario; reads: end jitter 3, no artifact, 100 per copy).** All scenarios: two surviving copies A (2% SNPs from the template) and B (3%) and the erased copy or copies; divergence of the erased copy from A about 2%. CTL: one erased copy E (the template), one isoform. ISO: E has a second isoform (internal exon 3 skipped), reads 67 : 33. SIB: two erased copies E1 (the template) and E2 (0.5% SNPs from it), both erased, 100 reads each (0.5% < the edge cutoff, so the frozen clustering joins them).
**Scoring (exact truth).** Expected transcripts: each (erased copy, isoform) with >= 5 reads in the net. A candidate RECOVERS a transcript T iff its core identity (end gaps excluded, edlib NW) to T is >= 0.999 and not lower than to the other expected transcripts of the family (a blend does not count). A candidate is SPURIOUS if no expected transcript reaches 0.999, REDUNDANT if it duplicates an already recovered one. Baseline = the frozen consensus per cluster.
**Bars (set now).** P1 CTL: PART splits at most 1 of the 12 control clusters and recovers every transcript the baseline recovers. P2 ISO: baseline recovers <= 60% of the expected transcripts (the minority isoform is lost, prediction), PART >= 80%. P3 SIB: baseline <= 60%, PART >= 80%. P4: spurious + redundant candidates <= 10% of PART's candidates over the three scenarios, and leaf purity (reads of a leaf that belong to its best-matching transcript) >= 95%.
**Real dev A, descriptive (not judged; H after, not blind):** clusters split, candidates added, and whether the best candidate against the erased copy raises the gate-aware identity x coverage >= 0.999 count (65 -> 80 with the metric alone) for the clusters with > 200 edits.
**Decision.** PART is adopted as an opt-in post-clustering step iff P1 to P4 hold; otherwise not adopted, with the failure.

#### Amendment 18 result (partition world seed 20261014; 36 families; 43 frozen clusters; 14 split; 57 candidates)

| scenario | expected transcripts | baseline (frozen consensus) recovers | PART recovers | clusters split / clusters | spurious + redundant candidates | leaf purity |
|---|---|---|---|---|---|---|
| CTL | 12 | 12 | 12 | 0 / 12 | 0 + 0 of 12 | 1.000 |
| ISO | 24 | 12 (50%) | **20 (83%)** | 9 / 12 | 0 + 1 of 21 | **0.849** |
| SIB | 24 | **19 (79%)** | **23 (96%)** | 5 / 19 | 0 + 1 of 24 | 0.958 |

**P1 met** (no control split, nothing lost). **P2 met** (baseline 50%, PART 83%). **P3 NOT met as written:** PART reaches 96% (bar >= 80%) but the baseline clause fails, baseline 79% not <= 60%, because the frozen clustering already separated 7 of the 12 sibling pairs. **P4 spurious + redundant met** (2 of 57 = 3.5%); **leaf purity NOT met for ISO** (0.849; bar 95%). **Decision: PART is not adopted by the registered rule.**
Cause of the ISO shortfall (3 of 12 clusters not split, all with a skipped exon of 271 to 322 bp; the 9 split include skipped exons of 93 to 310 bp): `minimap2 -x map-hifi` does not bridge a long deletion, it soft-clips the first part of the read (an iso1 read starts at consensus position 604 with a 286-base soft clip) so the skipped exon shows up as a change of coverage, not as a deletion run, and no variant column is found. This is an alignment limitation of the implementation, not of the partition logic.

---

### Amendment 19 (2026-10-08, after the Amendment 18 numbers, before any new number): global alignment for PART, fresh world

**Change.** PART aligns each read to the candidate with a semi-global edit-distance alignment (edlib `HW`: the read aligned in full inside the consensus, `task=path`, no distance cap), which keeps a long deletion as a deletion run, instead of `minimap2 -x map-hifi`. Reads are in transcript orientation (as the frozen pipeline's consensus is). Nothing else in PART changes (same variants, test, blocks, split rule, parameters).
**Test.** A fresh partition world, seed 20261015 (never run), E4 reads, 12 families per scenario. The bars of Amendment 18 are kept except that the baseline clause of P3 (a prediction about the frozen pipeline, not about PART) is replaced by PART >= baseline in every scenario: P1 (CTL: at most 1 of 12 clusters split, nothing lost); P2 (ISO: PART recovers >= 80%, and more than the baseline); P3 (SIB: PART >= 80%, and not fewer than the baseline); P4 (spurious + redundant <= 10% of PART's candidates over the three scenarios; leaf purity >= 95% in each scenario).
**Decision.** Adopted as an opt-in post-clustering step iff P1 to P4 hold on the fresh world. Real dev A is reported, not judged.

#### Amendment 19 result (fresh world seed 20261015; 47 frozen clusters, 13 split, 61 candidates; edlib alignment)

| scenario | expected | baseline recovers | PART recovers | split / clusters | spurious + redundant | purity |
|---|---|---|---|---|---|---|
| CTL | 12 | 12 | 12 | 0 / 12 | 0 + 0 | 1.000 |
| ISO | 24 | 12 | **21 (88%)** | 12 / 12 | 0 + **4** | **0.829** |
| SIB | 24 | 20 | **21 (88%)** | 1 / 23 | 0 + 3 | **0.875** |

P1, P2 and the recovery clauses of P3 met; **P4 NOT met** (spurious + redundant 7 of 61 = 11.5% against 10%; leaf purity 0.829 and 0.875 against 95%). **Decision: not adopted.** Diagnosis (unit test, not a result): a unit-cost edit-distance alignment prefers a messy alignment of 40 mismatches to a clean 100-base deletion, so the isoform shows up through noise correlated across the reads of one isoform, not as a deletion run; the split works by accident and assigns reads loosely, hence the purity.

---

### Amendment 20 (2026-10-08, after the Amendment 19 numbers, before any new number): spliced alignment, N counted as a deletion

**Change.** PART aligns with `minimap2 -ax splice:hq -uf --eqx` (checked on the failing cluster of Amendment 18: all 33 exon-skipping reads come out as ONE alignment with a > 100 bp N gap, where `map-hifi` clips) and counts an N gap as a deletion run in both the pileup and the read-by-variant table (for PART only; the Amendment 9 polish is unchanged). Nothing else changes. **Test:** a fresh world, seed 20261016, same bars as Amendment 19 (P1 to P4). **Decision:** adopted as an opt-in post-clustering step iff P1 to P4 hold.

#### Correction to the Amendment 18, 19 and 20 numbers (2026-10-08): a scoring bug in `parteval`, found by cross-checking the leaf composition

After the Amendment 20 run the leaf composition (`partition.json`) was perfect (all 12 ISO clusters split into a pure iso0 leaf of 67 reads and a pure iso1 leaf of 33) while the scorer reported purity 0.749 and 5 redundant candidates. The scorer took, for each candidate and expected transcript, the MAXIMUM core identity over the two orientations. A core identity is computed between the first and last run of 8 matches, so a garbage alignment (the wrong orientation, or an isoform against the other isoform) with one lucky run scores 1.0 over a tiny core. Fixed: the orientation is chosen by the GLOBAL identity and a core that covers less than half of the shorter sequence is `None` (tests added). Re-scoring the stored outputs of the three worlds (the partition output itself is unchanged):

| world (alignment) | CTL split / recovers | ISO baseline -> PART | SIB baseline -> PART | spurious + redundant of PART's candidates | purity CTL / ISO / SIB |
|---|---|---|---|---|---|
| 20261014 (Amendment 18, `map-hifi`) | 0 / 12 of 12 | 12 -> 21 of 24 (9 of 12 split) | 19 -> 24 | 0 of 57 | 1.00 / 0.905 / 1.00 |
| 20261015 (Amendment 19, edlib) | 0 / 12 of 12 | 12 -> 24 | 23 -> 24 | 1 of 61 | 1.00 / 0.994 / 1.00 |
| **20261016 (Amendment 20, spliced, the registered test)** | **0 / 12 of 12** | **12 -> 24 of 24 (12 of 12 split)** | 24 -> 24 | **0 of 60** | **1.00 / 1.00 / 1.00** |

**Verdicts with the corrected scorer.** Amendment 18 stays NOT adopted (P3 baseline clause: baseline 79% not <= 60%; ISO purity 0.905; `map-hifi` clips long deletions). Amendment 19 would have met P1 to P4 (the original "not adopted" came from the scoring bug; edlib's unit-cost alignment is still a poor aligner for long gaps, and its ISO split worked through correlated noise). **Amendment 20 meets P1 to P4 on the fresh world: PART is adopted as an opt-in post-clustering step** (`partition.py`, `PART_ALIGN=splice`). In the Amendment 20 world the frozen clustering had already separated every sibling pair (24 of 24 at baseline), so that world does not test the sibling gain; the 20261014 and 20261015 worlds, where the frozen edge rule joined 5 and 4 sibling pairs, show it (19 -> 24 and 23 -> 24 with 5 and 1 splits).

#### Amendment 18/20 real beds, descriptive (PART with `splice:hq -uf`, N as deletion; candidates aligned to the unmasked genome; the "best candidate" picks, per cluster, the candidate that matches the erased copy best, so it is an oracle upper bound, not a pipeline output)

| | dev A (212 clusters on the erased copy) | H (40; not blind) |
|---|---|---|
| clusters split | 72 (34%); leaves per split cluster 2: 33, 3: 9, 4: 8, 5: 9, 6: 6, 7: 5, 12: 1, 14: 1 | 19 (48%); up to 12 leaves |
| gate-aware identity x coverage >= 0.999: frozen consensus -> best candidate | 80 -> **106** | 17 -> **26** |
| the split clusters only | 40 -> 66 | 9 -> 18 |
| edit distance to the erased copy, split clusters | 579 -> 239 | 195 -> 10 |
| the clusters with > 200 edits (the original motivation) | 22 of 212, **split 0**, best candidate at 0.999: 0 | 3, split 0, 0 |

PART finds structure in a third to a half of the real clusters and, taking the best candidate, raises the count at the registered threshold by 26 (A) and 9 (H); the cluster-level evidence is a real mixture of members. It does NOT touch the 22 clusters with > 200 edits: their reads are not a mixture separable by linked variant columns (they are far from the assembled copy as a whole, probably reads that belong to another copy or isoform structure than the one the genome label points at). Unknown on real data: how many of the 417 candidates are real transcripts (no truth), and whether clusters with 5 to 14 leaves are over-split (the synthetic control splits 0 of 12, but real reads carry correlated errors and read-end heterogeneity the simulation lacks).

---

### Amendment 21 (2026-10-08, before the H number): the 22 dev-A clusters with more than 200 edits, and a truth-free flag for them

**What dev A shows (diagnosis, looked at before this amendment).** The 22 clusters are tiny (3 to 12 reads, median 5) and their reads are individually fine: baseline alignment to the unmasked genome has median divergence 0.0013 and MAPQ 60, as for all other clusters. What is broken is the consensus: the median read-to-consensus divergence (reads aligned to their own consensus, `splice:hq -uf`) is 0.016 (quartiles 0.009, 0.039) against 0.0015 for the other small clusters and 0.0014 for the large ones. The reads of such a cluster differ strongly in length and span (e.g. 667, 3,657 and 4,508 bases, genome starts up to 41 kb apart): the edge rule (divergence <= 0.00958 and a block covering half of the shorter read) lets a short read bridge two long reads that do not overlap each other, the connected component is a chain of partial or alternative-structure reads, and a global POA consensus of such a chain is a blend, not a transcript.
**Flag Q (truth-free).** A cluster's consensus is UNSUPPORTED if the median divergence of its reads aligned to it exceeds 0.00958 (the allele cutoff the edge rule already uses; no new constant). On dev A: 24 of 212 clusters are flagged, 16 of the 22 with > 200 edits and 8 others (sizes 3 to 75, edits 1 to 151); the 6 unflagged big-edit clusters have a consensus that explains its reads (divergence 0.001 to 0.009) but differs from the erased copy as a whole, which Q cannot see.
**Test on bed H (not used for the diagnosis).** Bars: Q flags >= 2 of the 3 H clusters with > 200 edits and at most 15% of the other H clusters (dev A: 4%). Reported: the flagged clusters' sizes.
**Decision.** Q is adopted as a quality flag on the cluster consensus ("consensus unsupported") iff both bars hold; it changes no consensus.

#### Amendment 21 result (2026-10-08; `consensus_support.py`)

| | dev A (design) | H (test) |
|---|---|---|
| clusters on the erased copy | 212 | 40 |
| flagged unsupported | 24 (sizes 3 to 75; 17 of them have 3 to 6 reads) | 3 (sizes 3, 8, 14) |
| of the clusters with > 200 edits flagged | 16 of 22 | **2 of 3** (bar >= 2) |
| other clusters flagged | 8 of 190 (4.2%) | **1 of 37 (2.7%)** (bar <= 15%) |

**Both bars met: Q is adopted as a quality flag ("consensus unsupported"); it changes no consensus.** What is explained: 16 of the 22 big-edit clusters on dev A (and 2 of 3 on H) are chains of 3 to 12 reads of very different spans whose POA consensus does not explain its own reads (median read divergence 0.016 against 0.0015). What is not: 6 big-edit clusters on dev A (and 1 on H) whose consensus explains its reads (divergence 0.001 to 0.009, 3 to 12 reads) but differs from the erased copy by 238 to 614 edits; their reads align to the unmasked genome individually at median divergence 0.0013, so the cause (reads from another copy or structure than the one the genome label points at, or a consensus that blends two structures the reads cannot expose) is not identified. In read terms the 22 clusters are 126 of about 59,000 reads (0.2%); in copy terms 10% of the clusters.

---

### Amendment 22 (2026-10-08, before any number from the new rule): chain-aware clustering, the edge must be a proper overlap

**What the pairs show (dev A, looked at before this amendment).** In the 16 flagged big-edit clusters, aligning the reads against each other pairwise (`splice:hq -uf`, the frozen edge rule applied to each pair) confirms only 104 of 216 pairs; 5 of the 16 clusters are not even connected by their own pairwise graph (3 have three isolated reads). The clusters were built from `map-hifi` read-to-seed edges, and the frozen edge accepts an alignment whose block covers only half of the shorter read (`block >= 0.5 x shorter`), wherever it sits: two reads that share a middle segment but continue differently at both ends are joined, and three such reads form a chain.
**Rule.** The frozen edge (gap-compressed divergence <= 0.00958 and block >= half of the shorter read) PLUS a proper overlap: at each end of the alignment at least one of the two reads must be reached, i.e. min(unaligned remainder of the query at that end, unaligned remainder of the target at that end) <= 0.00958 x block length (about 29 bases for a 3 kb alignment, the size of read-end ragged edges and of the 5' G run). A suffix-prefix overlap or a containment passes; an alignment that stops in the interior of both reads (a shared middle) does not. For a reverse-strand alignment the query start is paired with the target end. Nothing else changes (seed rounds, union-find components, minimum cluster of 3).
**Test.** Re-cluster dev A (design) and H (test) from scratch with the rule (`run_bed.py --proper --tag chain`), abPOA consensus, the Amendment 21 flag and the registered gate-aware metric on the clusters that overlap an erased copy. Bars, on H (test): (B1) reads of the erased copies clustered lose at most 3 points of coverage (frozen 95.7%); (B2) cluster purity by family stays 1.000; (B3) clusters with an unsupported consensus (Amendment 21) <= 1 of the 3 frozen; (B4) the gate-aware identity x coverage >= 0.999 count is not lower than the frozen 17, and the clusters with > 200 edits are <= 1 of 3; (B5) clusters on the erased copy do not fall by more than 5% (frozen 40). Dev A (design, direction only): unsupported <= 12 of 24, > 200 edits <= 11 of 22, fidelity count >= 80, coverage within 3 points.
**Decision.** Adopted as the cluster edge rule (opt-in flag `--proper`) iff B1 to B5 hold on H; otherwise not adopted with the failure.

#### Amendment 22 result (2026-10-08; `run_bed.py --proper --tag chain`, `eval_clusters.py`)

| | H frozen | **H chain (test)** | dev A frozen | dev A chain (design) |
|---|---|---|---|---|
| clusters | 47 | 54 | 221 | 312 |
| erased-copy reads clustered | 95.7% | **94.0%** | 96.4% | 95.6% |
| purity by family | 1.0000 | **1.0000** | 1.0000 | 1.0000 (1 background read in a family cluster) |
| clusters on the erased copy | 40 | **49** | 212 | 304 |
| unsupported consensus (Amendment 21) among them | 3 | **1** | 24 | 22 |
| gate-aware identity x coverage >= 0.999 | 17 (43%) | **26 (53%)** | 80 (38%) | 143 (47%) |
| clusters with > 200 edits | 3 | **2** | 22 | 27 |

B1 met (coverage -1.7 points), B2 met, B3 met (1 of 3), B5 met (+23%), **B4 NOT met: the fidelity count rises (17 to 26) but the clusters with > 200 edits are 2, bar <= 1.** Dev A direction bars fail too (unsupported 22 against <= 12; > 200 edits 27 against <= 11). **Decision: not adopted by the registered rule.** What the numbers say: the proper-overlap edge makes more, smaller clusters (+23% on H, +43% on A) and the share of them that reaches 0.999 on the erased copy rises about 10 points on both beds, but the chain survives where a short read is CONTAINED in two long reads that disagree with each other: each containment is a proper overlap, so A-B and B-C are edges although A and C are not compatible. A pairwise edge rule cannot see that; it needs the transitive check (Amendment 23).

---

### Amendment 23 (2026-10-08, after the Amendment 22 numbers, before any new number): star clustering inside each component

**Change.** Inside every component of the proper-overlap clustering with at most 60 reads: all-vs-all `minimap2 -x map-hifi -c`, a pair is COMPATIBLE iff it passes the Amendment 22 edge (divergence <= 0.00958, block >= half the shorter read, proper overlap); then LONGEST-FIRST STAR clustering: the longest unassigned read is the representative, its cluster is itself plus the unassigned reads compatible with it, repeat; clusters of fewer than 3 reads are dropped (their reads leave the clustering). A read joins only a cluster whose representative it is compatible with, never through another member, so a short read contained in two incompatible long reads goes to the longer one. Components with more than 60 reads keep the Amendment 22 result.
**Test (a, blind, exact truth): a synthetic CHAIN world, seed 20261017.** 24 families: 12 CHAIN (erased copy E with two equally expressed isoforms, iso0 full and iso1 skipping internal exon 4; 16 reads per copy, 40% of the reads 5'-truncated by up to 70%, so that fragments starting after exon 4 are contained in both isoforms) and 12 CTL (one isoform, same depth and truncation). Frozen clustering, Amendment 22 clustering and Amendment 23 clustering are scored the same way as Amendment 18: expected transcripts = (erased copy, isoform) with >= 5 reads in the net; a candidate recovers a transcript if its core identity (orientation by global identity) is >= 0.999 and not lower than to the other expected transcript; purity = reads of a cluster that belong to its best-matching transcript. **Bars:** CHAIN: Amendment 23 recovers >= 80% of the expected transcripts with purity >= 90%, and more than the frozen clustering; CTL: Amendment 23 recovers every transcript the frozen clustering recovers, and loses at most 5 points of read coverage; coverage over both scenarios >= 85%.
**(b) Real dev A and H, descriptive (both used for Amendment 22):** counts as in Amendment 22.
**Decision.** Adopted as the cluster rule (opt-in) iff the bars of (a) hold; (b) is reported.

#### Amendment 23 result (chain world seed 20261017; 12 CTL + 12 CHAIN families; frozen clustering, proper-overlap clustering, star clustering)

| | CHAIN expected / recovered | CHAIN reads clustered, purity | CTL expected / recovered | CTL reads clustered |
|---|---|---|---|---|
| frozen | 24 / 12 | 191 of 191, 0.50 | 12 / 12 | 192 of 192 |
| Amendment 22 (proper-overlap edge) | 24 / 12 | 191 of 191, 0.50 | 12 / 12 | 192 of 192 |
| **Amendment 23 (star step)** | 24 / **12 (50%)** | 182 of 191 (95%), **0.52** | 12 / 12 | 192 of 192 (100%) |

**Bar NOT met (CHAIN 50% against >= 80%, purity 0.52 against 90%): not adopted.** The CTL family is untouched (nothing lost, no coverage lost), the star step removes only 9 reads. Diagnosis: a full-length iso0 read and a full-length iso1 read (internal exon 4 skipped) align to each other in ONE alignment with a single 256 to 261 base deletion (block 2,122 of 2,124 bases), whose gap-compressed divergence `de` is 0.003: the pair passes the frozen rule, the proper-overlap test and the star compatibility, so the two isoforms stay together. The gap-compressed divergence hides exon-scale deletions by design.

---

### Amendment 24 (2026-10-08, after the Amendment 23 numbers, before any new number): compatibility counts every gap

**Change.** In the star step the pair's divergence is the full edit distance over the alignment block, NM / block (every inserted or deleted base counts), instead of the gap-compressed `de`; the pair is compatible iff NM / block <= 0.00958 (the same constant), the block covers half of the shorter read and the overlap is proper. A 260-base deletion in a 2,100-base block is 12%, far above the cutoff; HiFi indel errors (1 to 2 bases) are not. Nothing else changes (the seed rounds and the Amendment 22 edge keep `de`; the star step is applied to components of at most 60 reads).
**Test.** A fresh chain world, seed 20261018, same design and the same bars as Amendment 23: CHAIN: star clustering recovers >= 80% of the expected transcripts with purity >= 90% and more than the frozen clustering; CTL: recovers every transcript the frozen clustering recovers and loses at most 5 points of read coverage; coverage over both scenarios >= 85%. Real dev A and H reported, not judged (both used before).
**Decision.** Adopted as the cluster rule (opt-in) iff the bars hold.

#### Amendment 24 result (fresh chain world seed 20261018; real beds descriptive)

| chain world | CHAIN expected / recovered | CHAIN reads clustered, purity | CTL expected / recovered | CTL reads clustered |
|---|---|---|---|---|
| frozen | 24 / 12 | 192 of 192, 0.50 | 12 / 12 | 192 of 192 |
| **proper-overlap edge + star step (every gap counted)** | 24 / **24 (100%)** | 191 of 192 (99%), **0.91** | 12 / 12 | 192 of 192 (100%) |

**All bars met** (CHAIN 100% against >= 80% and more than the frozen 50%; purity 0.91 against >= 90%; CTL nothing lost; coverage 99 to 100% against >= 85%). **Adopted as the cluster rule, opt-in:** `run_bed.py --proper` then `run_chain.py` (star step, components of at most 60 reads). Purity 0.91 is the ambiguity of contained fragments, whose bases are identical in both isoforms and which go to the longer representative.
Real beds, descriptive (both used before; frozen -> proper edge -> + star step):

| | H frozen | H proper | **H + star** | dev A frozen | dev A proper | **dev A + star** |
|---|---|---|---|---|---|---|
| clusters on the erased copy | 40 | 49 | **54** | 212 | 304 | **298** |
| erased-copy reads clustered | 95.7% | 94.0% | **92.8%** | 96.4% | 95.6% | **94.9%** |
| purity by family | 1.0000 | 1.0000 | **1.0000** | 1.0000 | 1.0000 | **1.0000** |
| unsupported consensus (Amendment 21) | 3 | 1 | **0** | 24 | 22 | **8** |
| clusters with > 200 edits | 3 | 2 | **0** | 22 | 27 | **7** |
| gate-aware identity x coverage >= 0.999 | 17 (43%) | 26 (53%) | **29 (54%)** | 80 (38%) | 143 (47%) | **156 (52%)** |

The chain problem is largely resolved on the real beds: the clusters with more than 200 edits fall from 22 to 7 (dev A) and from 3 to 0 (H), the unsupported consensus sequences from 24 to 8 and from 3 to 0, and the share of erased-copy clusters at 0.999 rises 12 to 14 points; the cost is 1.5 to 2.9 points of read coverage (reads of unsupported chains leave the clustering) and +40% clusters. The 7 dev-A clusters still above 200 edits are the previously unexplained kind (consensus fits its reads, not the erased copy).

---

### Errata and limitations after the independent review (2026-10-08)

Four fresh reviewers (a claims-vs-data audit and three code/method reviews of Amendments 9 to 24) worked from a frozen read-only copy. What they found, what was fixed, and what was not. Result logs for the print-only analyses are in `docs/unmapped_rescue_results/`.
**Fixed in code (tests added):**
1. *Star step all-vs-all* (Amendment 23/24): `chain.py` ran minimap2 with `-N 20` and the default `-p 0.8`, which drops edges in components of 20 or more reads and fragmented ordinary clusters (a depth-50 control world: 31 clusters instead of 12, 92% of reads, registered control bar failed). Now `-N 500 -p 0.1` plus an integration test. The registered test (depth 16, components below 20 reads) was unaffected; on fresh worlds with the fix, depth 16, 22 and 50: isoform transcripts recovered 24 of 24 each, controls intact (100% of reads, 12 clusters), purity 0.89 to 0.93. Real beds re-run (`chain3`): H 58 clusters, 93.0% of erased-copy reads clustered, 55 on the erased copy, unsupported 0, clusters over 200 edits 0, gate-aware 0.999: 31; dev A 304, 95.0%, 294, unsupported 14, over 200 edits 13, 152. **The earlier dev-A figures (22 -> 7, 24 -> 8) came from the defective step; the corrected ones are 22 -> 13 and 24 -> 14.**
2. *Partition scorer* (`parteval`): the first fix still accepted fragments (a 3' half of a transcript scored 1.0 against it and, by tie order, against the other isoform). Now `ends.recovered_identity` requires the core to cover 90% of the transcript and of the candidate, ties are deterministic; the three partition worlds re-scored give the SAME tables (leaf composition had been checked independently), so verdicts stand.
3. *T4 search limit* (a templated run of 7 or more G was cut back to 6): the limit now follows the consensus run; test added.
4. *Unit tests*: `test_partition.py` built its fixture from a fresh `Random(7)` per character (a homopolymer); fixed. The edlib test asserted the opposite of the documented limitation (unit-cost edlib does not keep a 100-base deletion); it now documents it.
5. *Accounting of templated bases lost by T3* (Amendment 15): a clip that began with G but was not pure G counted as zero unaligned G. Corrected: at most 1 base lost in 157 of 197 trims (80%; H 28 of 37), still under the 90% bar; trims with no unaligned G at all are 61 (not 105). Verdict unchanged.
6. *Start-tolerant exactness* (Amendment 16): the modal-read-start convention scores a consensus that begins exactly at the true transcript as +1 when the jitter is one base (the untouched control is 21 of 24 exact). With a start-tolerant measure T4 is exact in 23 of 24 clusters, removes templated bases in 0, leaves an artifact G in 1 (this world and a fresh one, seed 20261301: 23, 0, 1; T3 11/4/9 and 15/1/8). T4 meets the synthetic bars (2) to (5); it remains not adopted because the dev-A bar fails (under-trim 36%; 23 of the 36 are clusters of fewer than 5 reads where T4 does nothing).
**Corrections to claims:**
7. *Amendment 17* was not tested on data "not used for this rule": the worlds had been scored for T3/T4, whose trim is the same operation, and the guard against false passes cannot fire on error-free consensus sequences. Adopted as a reporting CONVENTION; it changes the registered O3 flag (LRPAP1 c01 and the c00 allele cluster flip to flagged). A truly fresh world (seed 20261301) and the reviewer's (20261101) agree: 19 and 14 -> 24 of 24, controls identical, no false passes.
8. *Amendment 24 purity bar* (0.90): across eight chain worlds the purity is 0.85 to 0.93 (below 0.90 in four with the defective step; 0.89 to 0.93 with the fixed one). The rule is judged on the recovery and control bars; purity is limited by fragments that are identical in both isoforms. "All bars met" holds for the registered seed only.
9. *Remaining over-200-edit clusters* (Amendment 21/24): 71 of the 126 reads of the 22 clusters leave the clustering (dev A), 14 of 14 on H, so part of the improvement is removal; some remaining ones are unsupported (not "consensus fits reads"). Six of the 22 dev-A clusters hold reads of BOTH orientations (the pool graph ignores strand and the consensus is not oriented), an undisclosed mechanism for "consensus unsupported". The star step refines only components of at most 60 reads (4 to 8% of the reads); the mixed-isoform consensus in large clusters is untouched.
10. *Partition (Amendments 18 to 20)*: (a) validated only against independent errors; with errors correlated across a read subset (a hotspot of fixed columns in a quarter of the reads) PART split 10 of 12 control clusters in the reviewer's test, and two linked heterozygous SNPs split 5 of 12, so the real-bed split counts (34 to 48%, up to 14 leaves) are not evidence of distinct transcripts; (b) it is blind to an extra exon carried by a minority (an insertion is one variant, never a block) and its sensitivity floor depends on the variant length (a 280-base skip in 8 of 100 reads is missed); (c) the sibling gain (19 -> 24, 23 -> 24) comes from worlds that used the rejected alignments, and the adopted configuration separated 0.5% pairs in a forced-join test but only 10 of 12 at 0.2%; (d) the real-bed "80 -> 106" counts `best or frozen` (monotone by construction), 8 of the 26 gains are a shorter leaf scoring 1.0, and 28 of the 72 split dev-A clusters have two largest leaves of equal length differing by 1 to 15 substitutions (an allele signature in a diploid animal). "A real mixture of members" is an interpretation. PART stays opt-in and is not validated for correlated errors or alleles.
11. *Stale or wrong statements fixed in place:* "held-out H" (Amendment 2 already declared H not held out for the arms that followed; "not tuned on H" is what holds); "30 of 30 families" and "0 of 30" (26 clusters exist); "joined 5 and 4 pairs" (5 and 1); in Amendment 16 the j_hat tallies (6, 6, 6, 7 = 25) were not per erased copy; in Amendment 14 both failures were `too_long` (the `not_g` cluster passes); in Amendment 15 the design-check paragraph (39,351; 28,581) differs from the result file (39,403; 28,614) and "p = 1" for controls is the n < 20 sentinel; in Amendment 21 "17 of them" is 16 and the median divergence is 0.0176 (0.0094, 0.0419); in Amendments 18/20 "417 candidates" covers all 221 clusters (407 for the 212); in Amendment 5 GWFAM163 attracts six wrong clusters (125 wrong reads in 10 clusters) and two clusters are lost, the +243 holds; in section 4 "26 of 38" is 26 of 39. The real beds are real reads with SYNTHETIC erasure, and the library gate and artifact distribution are computed from H's own survivor reads.
**Known statistical looseness (no verdict changes):** the binomial gate tests p = 0.5, not 0.25; the polish test divides alpha by 3L but tests more alternatives per column (anti-conservative, so Amendment 9's refutation is unaffected); the artifact length distribution is measured from soft clips and undercounts artifact G's that match the genomic flank (a(0) 0.21 against 0.05 simulated), which biases the templated-count estimate upward.
**Not checkable from the frozen copy, now backed by logs or still open:** the real-bed trimming and metric analyses (logs saved), the margin sweeps on dev A, the O3 maternal study (no result files in the copy), Amendment 8's score-rule per-copy breakdown (stored for the divergence rule only).

---

### Amendment 25 (2026-10-08, before any number): specificity controls for PART and for the O3 flag

**Why.** The independent review showed that the synthetic controls so far have independent errors, haploid copies and a complete net, so they cannot expose the two ways the pipeline can call something that is not there: (i) PART splitting a cluster whose reads share a recurrent sequencing error (a HiFi homopolymer hotspot), and (ii) the O3 flag calling an allele of a present copy, or a family with nothing missing, an absent copy.
**Part A: a guard for PART (block purity).** A linked block is a real haplotype or isoform when reads carry (almost) all of its variants or (almost) none; a hotspot is carried partially and ALSO appears at a few percent in the reads that do not carry it. For a candidate block with m variants over n reads (carriers = at least half of the variants they cover): `Dev` = the number of (read, variant) pairs in which the read disagrees with its own group (a carrier lacking the variant, a non-carrier having it); `e_bg` = the non-consensus allele rate outside the variant columns (floor 1e-4); the block is REJECTED iff P(Binomial(number of pairs, e_bg) >= Dev) < 0.05 / (number of candidate blocks); the next block is tried. No other change to PART.
Bars: read-level stress harness (`stress_part.py`, real `splice:hq` alignment, 12 transcripts per setting, fresh seeds). (A1) Hotspot control (75 good reads with 3% error at F fixed columns, 25 low-quality reads with probability p): splits <= 1 of 12 for (F, p) = (20, 0.3), (20, 0.5), (40, 0.5); independent-error control 0 of 12. (A2) True structure retained: an exon skip in 33% of the reads, split in >= 11 of 12; two siblings 0.5% apart (6 or more differing columns), split in >= 11 of 12; two linked heterozygous SNPs, split (allelic heterogeneity is reported by PART, classified by Part B). (A3) The registered partition world is not harmed: a fresh world, seed 20261401 (`--partition-world`): CTL splits 0 of 12, ISO recovers >= 22 of 24.
**Part B: separating an allele from an absent copy in the O3 flag.** For a cluster consensus: R = identity x coverage (gate-aware) of its best hit on the reference genome, T = the same on the truth genome (the other haplotype assembly, or in a synthetic world the genome with every copy). Registered flag (unchanged): T >= 0.999 and R < 0.999. New classification: COPY = T >= 0.999 and R < 1 - 0.00958 (more than the allele cutoff from every reference locus); ALLELE = T >= 0.999 and 1 - 0.00958 <= R < 0.999; PRESENT = R >= 0.999. Only COPY is an O3 call.
Bars: synthetic specificity world (seed 20261402, `--spec-world`): 12 NULL families (three copies, all present), 30 ALLELE families (A and B present; the individual is heterozygous at A, the second allele A' differs from A by 0.25%, 0.5%, 1%, 2%, 4%, six families each, A' absent from the reference), 12 ERASED families (copy E at 2% from its nearest survivor, the positive control). (B1) NULL: no cluster is classified COPY and none is flagged by the registered rule. (B2) ALLELE: for A' at 0.25% and 0.5%, 0 COPY calls (the reads either do not enter the net or form clusters classified ALLELE); the table for 1%, 2% and 4% is REPORTED, not judged (an allele more than 1% from its other haplotype is indistinguishable by sequence from a paralog; it is what RNA alone cannot decide). (B3) ERASED: >= 95% of the erased-copy clusters are COPY. (B4) LRPAP1, a sanity check on known data: the c00 allele cluster is ALLELE, the c01 clusters are COPY.
**Decision.** The guard is adopted into PART iff A1 to A3 hold; the three-way classification is adopted for the O3 call iff B1 to B4 hold. Failures are reported as such.

#### Amendment 25 result (2026-10-08; `stress_part.py`, `flagmetric.o3_class`, `run_synth.py o3class`; logs in `docs/unmapped_rescue_results/`)

**Part A, the PART block-purity guard (read-level stress, 12 transcripts per setting, real `splice:hq` alignment).**

| setting | guard OFF (seed 5001) | **guard ON, seed 5001** | **guard ON, seed 5002** |
|---|---|---|---|
| independent-error control | 0 of 12 split | **0** | **0** |
| hotspot, F = 20 columns, p = 0.3 | 10 of 12 | **1** | **0** |
| hotspot, F = 20, p = 0.5 | 12 of 12 | **0** | **0** |
| hotspot, F = 40, p = 0.5 | 12 of 12 | **0** | **0** |
| exon skip in 33% of the reads | 12 of 12 | **12** | **12** |
| two siblings 0.5% apart | 12 of 12 | **12** | **12** |
| two linked heterozygous SNPs | 6 of 12 | 6 | 11 |

A1 met (hotspots <= 1 of 12 in every setting, control 0). A2 met for the skip and the siblings (12 of 12; the bar was >= 11); the heterozygous case is not a numeric bar and PART finds two linked SNPs in 6 and 11 of 12 transcripts (a 2-column block is at the edge of its power). A3 met on a fresh partition world (seed 20261401): CTL 0 of 12 split, ISO 24 of 24 recovered (12 of 12 split), SIB 24 of 24, purity 1.00, no spurious or redundant candidate. **The guard is adopted into PART.** Effect on the real beds (descriptive, same reads): dev A 63 of 221 clusters split (73 before), 330 candidates (417), at most 7 leaves (14), best-candidate gate-aware 0.999 count 103 (106); H 17 of 47 (20), 88 candidates (118), at most 10 leaves (12), 24 (26). The guard removes the many-leaf splits; the remaining splits are still not validated against alleles (see the Errata). Reproducibility note: two H clusters' identical consensus sequences gave alignments that differ by one edit between two minimap2 batches, which moves the H frozen count 17 -> 16.
**Part B, the allele-vs-copy classification (specificity world seed 20261402; E4 reads; library gate closed in this world).**

| scenario | true kind | reads in net | clusters | COPY | ALLELE | registered flag |
|---|---|---|---|---|---|---|
| NULL (3 copies, nothing missing) | present | 0 | 0 | 0 | 0 | 0 |
| ALLELE 0.25% | allele | 0 / 480 | 0 | 0 | 0 | 0 |
| ALLELE 0.5% | allele | 7 / 480 | 1 | **0** | 1 | 1 |
| ALLELE 1% | allele | 408 / 480 | 6 | 4 | 2 | 6 |
| ALLELE 2% | allele | 480 / 480 | 6 | 6 | 0 | 6 |
| ALLELE 4% | allele | 480 / 480 | 6 | 6 | 0 | 6 |
| ERASED (positive control) | copy | 960 / 960 | 12 | **12** | 0 | 12 |

B1 met (no cluster at all from the 12 NULL families, so nothing is flagged), B2 met (0 COPY calls for alleles at 0.25% and 0.5%; at 0.5% the one cluster is ALLELE where the registered flag calls it a copy), B3 met (12 of 12 erased copies are COPY). **B4 met on LRPAP1:** the 21-read c00 cluster (mat best 0.9959) is ALLELE, the 283-read and 4-read c01 clusters (mat best 0.9899, 0.9895) are COPY; the 3-read cluster has no call. Reported, not judged: an allele 1% from its other haplotype is already called COPY in 4 of 6 families and 2% or more in all (6 of 6): by sequence alone RNA cannot tell an allele that far apart from a paralog, which is the limit of any RNA-only O3 call; the c01 consensus sits only 0.0005 below the allele cutoff on LRPAP1 (0.9899 against 0.99042). **The three-way classification is adopted for the O3 call; the registered flag (COPY or ALLELE) is kept as the screening statistic.**

---

### Amendment 26 (2026-10-08, before any classification): the three-way O3 call on real data with an independent truth

**Setting.** The maternal-reference study (`O3_MATERNAL_REFERENCE_2026-10-08.md`) has real KB3781 reads, the animal's own two haplotype assemblies and, for each judged copy, the identity of its nearest relative on the reference haplotype (section 2 tables). The reference for the call is one haplotype (R), the truth the other (T). Pool = every one of the 35,094 reads of the 34 families and LRPAP1 whose primary alignment on the reference haplotype is unmapped or has divergence > 0.00958; clusters by the adopted rule (`--proper` edge, star step with the fixed all-vs-all, minimum 3 reads); abPOA consensus; R and T = gate-aware identity x coverage of the consensus' best hit on each haplotype (`flagmetric.o3_class`, library gate computed from this dataset's own shared reads).
**Prediction, written from the study's tables before any consensus is classified.** With the MOTHER as reference the only judged copy with no relative within the allele cutoff is GWFAM175_B0 (nearest relative 0.9126); every other judged copy has one at 0.9911 to 1.0000 (GWFAM208_B0 1.0000, 214_B0 0.9995, 205_B0 0.9935, 64_B17 0.9999, 64_B11 0.9980, 227_B0 0.9911, LRPAP1_p12 0.9959). So COPY is expected for GWFAM175_B0 and, from the LRPAP1 example, for the LRPAP1 c01 cluster; ALLELE or PRESENT for the others. With the FATHER as reference all 13 judged copies have a relative at 0.9927 or closer (GWFAM382_B0 matches neither haplotype within 5%): no COPY call is expected.
**Bars (set now).** (R1) Mother as reference: at least one cluster of >= 20 reads whose majority label is GWFAM175_B0 is COPY. (R2) Mother as reference: no cluster whose majority label is one of GWFAM208_B0, 214_B0, 205_B0, 64_B17, 64_B11, 227_B0 or LRPAP1_p12 is COPY, and no cluster of 'shared' reads (reads that map to the same place on both haplotypes), LRPAP1 reads excluded, is COPY. (R3) Father as reference: zero COPY calls among the clusters of the judged-copy labels, and at most 1 COPY call in total. (R4) LRPAP1 in the same flow (mother as reference): the c01 clusters are COPY and the c00 cluster is ALLELE. Reported, not judged: every COPY call in each direction with its reads, labels and identities; the calls the registered flag would make; the number of ALLELE and PRESENT clusters; (for bed H, descriptive) the three-way classes of the erased-copy clusters against the masked reference.
**Decision.** The three-way call is reported as supported on real data iff R1 to R4 hold; each failure is reported as such.

#### Amendment 26 result (2026-10-08; `o3_haplotype.py`, logs `docs/unmapped_rescue_results/o3_haplotype.{mat,pat}.eval.log`)

Pools: mother as reference 2,888 net reads (GWFAM175_B0 283 of its 304 reads, 'shared' 2,392, LRPAP1 included), 37 clusters of >= 3 reads; father as reference 2,824 net reads, 36 clusters. The library gate is open in both (G 20,247 against A, C, T 34 pure 5' clips).

| mother as reference, truth father | clusters | COPY | ALLELE | PRESENT | no call | registered flag |
|---|---|---|---|---|---|---|
| GWFAM175_B0 | 1 (269 reads) | **1** | 0 | 0 | 0 | 1 |
| LRPAP1 | 6 | 2 (the 282-read and the 4-read c01 clusters) | 0 | 0 | 4 | 2 |
| 'shared' reads | 28 | **5** (356, 63, 28, 8, 65 reads) | 6 | 0 | 17 | 11 |
| 'ambiguous' | 2 | 0 | 0 | 0 | 2 | 0 |

Father as reference: 36 clusters; COPY 10 (7 'shared' with 210, 79, 134, 480, 3, 4 and 101 reads; 2 'ambiguous'; 1 LRPAP1); the judged-copy labels GWFAM70_B4 ALLELE, GWFAM26_B2 and GWFAM382_B0 no call.
**Bars. R1 met** (the GWFAM175_B0 cluster is COPY, R = 0.9306). **R2 NOT met** (no judged-locus cluster is COPY, but 5 of the 28 'shared' clusters are). **R3 NOT met** (no COPY among the judged-copy labels, but 10 COPY calls in total against <= 1). **R4 partly:** both c01 clusters are COPY, but the c00 allele cluster of the earlier flow is not formed in this pool (the four small c00 clusters have T below 0.999: no call), so ALLELE is not shown. **The three-way call is NOT supported on real data by the registered bars.**
What the failures are: every one of the 'shared' COPY clusters has its reference hit and its truth hit on the SAME chromosome of the two haplotypes (chr5_mat 49.05 Mb with chr5_pat 41.47 Mb, chr8, chr11, chr2, chr1), i.e. it is a diverged counterpart at the corresponding locus: R = 0.973 to 0.990 (and 0.83 to 0.93 in the father's direction), T = 1.0. Real haplotypes of this animal differ by 1 to 7% at a minority of expressed loci (5 of 28 'shared' clusters with the mother as reference, 7 of 24 with the father), above the 0.958% allele cutoff, so a purely sequence-based call over-calls. The GWFAM175_B0 hit (chr5_mat 47.66 Mb) sits at the same offset to its pat locus (7.64 Mb) as the neighbouring 'shared' cluster (7.58 Mb): its nearest mother's relative is at the orthologous position too, 7% away. What separates the two kinds is not distance but whether the T locus is the reference locus' best match in return (one-to-one orthology), Amendment 27.

---

### Amendment 27 (2026-10-08, after the Amendment 26 numbers, before any reciprocal result): a COPY must have no reciprocal best match

**Rule.** For a COPY-called cluster, take the reference locus X of its best hit and build the X TRANSCRIPT, the reference segments the consensus alignment covers (exons joined, introns removed, from the `cg` CIGAR of the existing PAF). Map the X transcript to the TRUTH haplotype (`splice:hq -uf -N 5`, the same index). X is the ortholog of the truth locus Y0 iff the best record of the X transcript (by matches) lies on the same chromosome and overlaps Y0 (the truth hit of the consensus) by >= 50% of the shorter. Reciprocal COPY -> ALLELE ("diverged ortholog"); non-reciprocal -> COPY confirmed (the reference locus belongs to another truth locus, so Y0 has no counterpart). No new constant.
**Predictions (from the study's own tables, before the test).** Mother as reference: the GWFAM175_B0 cluster and the LRPAP1 c01 clusters are non-reciprocal (the study found that their nearest relatives are claimed by other copies); the 5 'shared' COPY clusters are reciprocal. Father as reference: the 7 'shared' COPY clusters are reciprocal.
**Bars.** (S1) mother as reference: GWFAM175_B0 and the 282-read LRPAP1 c01 cluster stay COPY. (S2) at least 4 of the 5 'shared' COPY clusters (mother) and at least 6 of the 7 (father) become ALLELE. (S3) after the rule, COPY calls from 'shared' clusters are <= 1 in each direction. Reported: every call, the Y1 locus found for each X transcript, the 'ambiguous' and LRPAP1 COPY clusters in the father's direction.
**Decision.** The reciprocal rule is adopted into the three-way call iff S1 to S3 hold.

#### Amendment 27 result (2026-10-08; `reciprocal.py`, `o3_haplotype.py reciprocal`, logs `o3_haplotype.{mat,pat}.reciprocal.log`; `o3_bed.py`)

**Mother as reference** (8 COPY clusters): the GWFAM175_B0 cluster (269 reads) stays COPY (its mother's locus, chr5_mat 47.66 Mb, back-maps to a different father's locus, chr5_pat 39.97 Mb, at 0.932); the 282-read and the 4-read LRPAP1 c01 clusters stay COPY (the mother's locus is c00, whose father's counterpart is c00 on chr3, not c01 on chr12); **all 5 'shared' COPY clusters are reciprocal** (356, 63, 28, 8, 65 reads: each back-maps onto its own father's locus at 0.977 to 0.989) and become ALLELE (diverged orthologs). S1 met, S2 met (5 of 5), S3 met (0 'shared' COPY).
**Father as reference** (10 COPY clusters): 4 of the 7 'shared' COPY clusters are reciprocal (210, 79, 134, 101 reads) and become ALLELE, as does the 5-read LRPAP1 cluster; **3 'shared' clusters stay COPY** (the 480-read cluster, whose father's locus chr5_pat 39.97 Mb back-maps to a DIFFERENT mother's locus, chr5_mat 47.55 Mb, at identity 1.0000; two chr1 clusters of 3 and 4 reads) and the two tiny 'ambiguous' clusters stay COPY. **S2 NOT met (4 of 7 against >= 6), S3 NOT met (3 against <= 1). The reciprocal rule is NOT adopted by the registered bars.**
What the residual calls are: the 480-read cluster is the mirror image of the mother-direction GWFAM175_B0 call. In the mother's direction the father's copy at 40.03 Mb has no ortholog (its nearest mother's locus, 47.66 Mb, 93% identical, is claimed in return by a different father's locus); in the father's direction the mother's copy at 47.66 Mb has no ortholog (its nearest father's locus is 93% identical and claimed by the father's other copy). The two haplotypes of this animal carry two paralogs of the GWFAM175 tandem array that the other lacks: a real copy-number polymorphism, expressed in 269 and 480 reads. The study's synteny-based truth list (loci absent from the father/mother) lists the first and labels the second 'shared'; the rule found a true absent copy the registered truth misses, not a false call. The bars assumed the list was complete (the study had already said it is not). The two chr1 clusters (3 and 4 reads, R = 0.83) and the two 'ambiguous' clusters are too small to judge.
**Bed H (erasure, masked genome as reference; descriptive):** 55 clusters on an erased copy: 31 COPY (the same 31 that reach 0.999), 24 no call (not reconstructed at 0.999), 0 ALLELE; 3 background clusters, no call. By construction the net holds only reads more than the allele cutoff from the reference, so an erased copy whose nearest survivor is within the cutoff leaves no cluster to call.
**Verdicts.** The three-way sequence-only call (Amendment 25) is supported where truth is exact (synthetic) and on bed H, and over-calls on real haplotypes where 5 of 28 (mother) and 7 of 24 (father) expressed 'shared' clusters are diverged orthologs 1 to 7% apart. The reciprocal rule removes 9 of those 12 and leaves one genuine copy-number polymorphism and four tiny clusters; by the registered bars it is not adopted, and the result is best read as: the call needs one-to-one orthology (not a distance) to separate alleles from copies, and real haplotypes show both a high rate of diverged orthologs and copies the synteny truth list does not contain.

---

### Amendment 28 (2026-10-08, before any result): discovery of reference-absent sequence from the transcriptome, the 959 reads unmapped on the primary assembly

**Question.** The aim is sequence the reference genome does not contain but the transcriptome reaches, with no truth built in. Real case: KB3781 Iso-Seq, the 959 reads that do not align to the combined primary assembly (`o3_mat/reads/R_unm.fa`). Independent confirmation is used only to grade the result: the animal's maternal and paternal assemblies (132 of the 959 reads map on the mother's, 3 on the father's, 125 of them in one 71-kb stretch of mat chromosome 12; 824 map on neither).
**Pipeline (no truth, no haplotype).** Cluster the 959 reads with the adopted rule (`--proper` edge, star step with the fixed all-vs-all, minimum 3 reads), abPOA consensus per cluster (<= 100 longest reads), and align each consensus (`splice:hq -uf -N 5`, gate-aware identity x coverage, the library gate computed from this dataset's shared reads, open) to the primary (R), to the mother's (Tm) and to the father's (Tp) assembly. Classes: CONFIRMED = R below 1 - 0.00958 and max(Tm, Tp) >= 0.999 (absent from the primary, present in the animal's other assembly); NOVEL = R below 1 - 0.00958 and max(Tm, Tp) < 0.999 (a reconstructed sequence none of the three assemblies holds); ALLELE-LIKE = R in [1 - 0.00958, 0.999); PRESENT = R >= 0.999.
**Control (specificity).** 959 reads sampled (seed 1) from `o3_mat/reads/R_LRP.fa` (LRPAP1 reads; every one of them maps on the primary), same pipeline, same classes.
**Predictions.** The 125-read stretch forms one cluster that is CONFIRMED against the mother's assembly; most of the other 824 reads stay unclustered; a few small clusters are NOVEL.
**Bars (set now).** (D1) the largest cluster (>= 100 reads) is CONFIRMED. (D2) at least 80% of the 135 reads that map on the mother's or father's assembly lie in CONFIRMED clusters. (D3) control: no CONFIRMED cluster and at least 98% of its clusters PRESENT or ALLELE-LIKE (no NOVEL). Reported, not judged: the number, reads and lengths of the NOVEL clusters, the reads left unclustered, the clusters of 'unmapped' reads that are PRESENT or ALLELE-LIKE.
**Decision.** The discovery run is reported as working iff D1 to D3 hold; every NOVEL cluster is a candidate for DNA verification (consensus k-mers in the individual's genomic reads), not a call.

#### Amendment 28 result (2026-10-08; `discover.py`; logs `docs/unmapped_rescue_results/discover.*.log`)

959 unmapped reads -> 3 clusters of >= 3 reads (113, 4 and 3 reads; 120 reads clustered, 839 unclustered). The 113-read consensus (5,694 bp) aligns to the MOTHER's assembly at chr12_mat 95.96-96.03 Mb with identity 1.0000 over 5,557 of 5,694 bases (no hit on the primary), gate-aware identity x coverage Tm = 0.9759: **CLASS NOVEL, so D1 NOT met and D2 NOT met (0 of 134 confirmable reads in CONFIRMED clusters; 113 are in a cluster).** Control (959 clean LRPAP1 reads): 14 clusters, 7 PRESENT (362 reads), 6 ALLELE-LIKE (34), **1 NOVEL (462 reads, 4,244 bp, R = 0.77)**, 0 CONFIRMED: **D3 NOT met** (13 of 14 = 93%).
Why. (1) The 137 unaligned 5' bases of the 113-read consensus are real sequence of the same locus, not an artifact: normal composition, and mapped on their own they fall in the mother's locus in two pieces (46 bases at 96,028,629 and 90 bases at 95,981,351): small terminal exons the spliced aligner did not join, so the registered best-record coverage undercounts. (2) The control NOVEL cluster is a 462-read component (above the 60-read limit of the star step) whose consensus is a blend of reads that overlap in different parts of a multi-copy locus; the Amendment 21 support test would reject it. So the classes missed their own aim through two defects of the definition, not of the pipeline. Registered as written, D1 to D3 fail.

---

### Amendment 29 (2026-10-08, after the Amendment 28 numbers, before any revised classification): terminal-exon rescue and a support requirement

**Changes to the discovery classification only (clustering, consensus and alignments are unchanged).** (a) TERMINAL-EXON RESCUE: the consensus bases left unaligned at either end of the best record (>= 20 bases) are mapped (`minimap2 -x sr -N 5`) to the genomic window [record start - 5 kb, record end + 5 kb] of the same assembly; pieces with identity >= 0.98 are added to the coverage (union of query intervals) and to the identity (matches over block); applied to all three assemblies; the 5' G rule of Amendment 17 is unchanged. (b) A cluster whose consensus is UNSUPPORTED (median divergence of its reads to it above 0.00958, Amendment 21) is class UNSUPPORTED and is never CONFIRMED or NOVEL.
**Predictions (from the diagnosis above).** The 113-read cluster becomes CONFIRMED (the two pieces account for the 137 bases); the control 462-read cluster becomes UNSUPPORTED; the 4- and 3-read clusters stay NOVEL or UNSUPPORTED.
**Bars.** D1: the largest cluster of the 959 reads is CONFIRMED. D2: >= 80% of the 134 confirmable reads are in CONFIRMED clusters (113 of 134 = 84% if D1 holds). D3: on the original control (seed 1, not blind now) and on a FRESH control (959 reads, seed 2, blind): no CONFIRMED cluster and no NOVEL cluster. Reported: the NOVEL and UNSUPPORTED clusters of the 959 reads with sizes and lengths.
**Decision.** The discovery run is reported as working iff D1 to D3 hold; every NOVEL cluster remains a candidate for DNA verification, not a call.

#### Amendment 29 result (2026-10-08; `discover.py`, `rescue.py`; logs `discover.unm.log`, `discover.control.log`)

**The 959 unmapped reads, revised classes:** 3 clusters (120 reads clustered, 839 unclustered). **The 113-read cluster (5,694 bp) is CONFIRMED:** no hit on the primary, Tm = 1.0000 after the terminal-exon rescue (the 137 unaligned bases are two small exons of the same locus on the mother's chr12), median read divergence to the consensus 0.0010 (supported). **D1 met. D2 met: 113 of the 134 confirmable reads (84%) are in the CONFIRMED cluster; the other 21 are in no cluster.** The other two clusters (4 reads, 2,445 bp; 3 reads, 2,330 bp) are supported (read divergence 0.0008 and 0.0009) and **NOVEL: no alignment at all on the primary, the mother's or the father's assembly**, 7 reads in total, candidates for DNA verification.
**Control (original, 959 LRPAP1 reads, seed 1): D3 NOT met.** 14 clusters: 7 PRESENT (362 reads), 6 ALLELE-LIKE (34), and the 462-read cluster, now supported (median read divergence 0.0013), is CONFIRMED with R = 0.964 on the primary and Tm = 0.9993 on the mother's: its consensus is the mother's haplotype of a 4.2 kb stretch where the combined primary carries the father's, 3.6% away. So the control was not clean: "every read maps on the primary" is not "every read matches the primary closely", and the animal's two haplotypes differ there. The pipeline did what it should (a sequence of the animal absent from the primary, present in an assembly of the same animal); the control design was wrong. Registered as written, D3 fails; a clean control is Amendment 30.

---

### Amendment 30 (2026-10-08, after the Amendment 29 numbers, before the clean control): a control that is clean on the reference

**Control.** 959 reads (seed 3) sampled from the 34-family set that carry the label 'shared' and whose primary alignment has divergence <= 0.00958 on BOTH the mother's and the father's assembly (a read clean on both haplotypes is clean on the mosaic primary), same pipeline, same classes. **Bar (D3').** No cluster is CONFIRMED, NOVEL or UNSUPPORTED, and at least 95% of the clustered reads lie in PRESENT or ALLELE-LIKE clusters. **Decision.** The discovery run is reported as specific iff D3' holds.

#### Amendment 30 result (2026-10-08; `discover.py control_clean`; log `discover.control_clean.log`)

Clean control (959 reads clean on both haplotypes, seed 3): 46 clusters, 660 reads clustered; PRESENT 31 (444 reads), ALLELE-LIKE 12 (195), **NOVEL 3 (21 reads: 14, 4, 3)**, CONFIRMED 0, UNSUPPORTED 0. 96.8% of the clustered reads are PRESENT or ALLELE-LIKE (bar >= 95% met) but **D3' NOT met** (no NOVEL cluster). The three NOVEL clusters are not absent sequences: every one has a full-length hit on all three assemblies, at identity x coverage 0.977 (2,729 bp, the first 50 bases unaligned and not rescued), 0.988 on the primary and 0.9990 on the father's (1,977 bp, 0.0001 under the bar) and 0.962 (1,117 bp, first 43 bases unaligned). So "NOVEL" lumped a sequence no assembly holds with a sequence that has a good but imperfect hit (a 0.001 cliff at the bar; 40 to 50 non-G bases at the 5' end). The 5' non-G end of 43 to 50 bases is a second kind of end artifact, not examined here.

---

### Amendment 31 (2026-10-08, after the Amendment 30 numbers, before the fresh control): tiers, as in the candidate list the user asked for

**Change.** The class that means "no assembly holds this" must say so: NOVEL = the best identity x coverage over the primary, the mother's and the father's assembly is below 0.9 (or no hit). A supported consensus with R below 1 - 0.00958 and a best truth-assembly score in [0.9, 0.999) is DIVERGED (a hit exists, nowhere at the bar: a paralog, a structural difference, or incomplete ends); CONFIRMED, ALLELE-LIKE, PRESENT and UNSUPPORTED as before. 0.9 is the boundary used for the candidate tiers (below 0.9 absent, 0.9 to 0.99 diverged, above 0.99 allele-like).
**Test.** The class change is applied to the three runs already made (the 959 unmapped reads, the LRPAP1 control, the clean control) and to a FRESH clean control (959 reads, seed 4, blind). **Bars:** the 959 reads: the 113-read cluster stays CONFIRMED and the two small clusters stay NOVEL; clean controls (seed 3 and seed 4): no CONFIRMED, NOVEL or UNSUPPORTED cluster and >= 95% of the clustered reads PRESENT or ALLELE-LIKE (the DIVERGED share is reported). **Decision.** The discovery run is reported as specific iff the control bars hold on the fresh control.

#### Amendment 31 result (2026-10-08; `discover.py`, logs `discover.{unm,control,control_clean,control_clean4}.log`; sequences `discovery_candidates.fa`)

| run | clusters / reads clustered | CONFIRMED | NOVEL | DIVERGED | ALLELE-LIKE | PRESENT |
|---|---|---|---|---|---|---|
| **959 unmapped reads** | 3 / 120 | **1 (113 reads)** | **2 (4 and 3 reads)** | 0 | 0 | 0 |
| LRPAP1 control (seed 1, not clean) | 14 / 858 | 1 (462) | 0 | 0 | 6 (34) | 7 (362) |
| clean control (seed 3) | 46 / 660 | 0 | 0 | 3 (21) | 12 (195) | 31 (444) |
| **fresh clean control (seed 4, blind)** | 54 / 684 | **1 (9 reads)** | 0 | 0 | 20 (351) | 33 (324) |

**The 959 reads: D1 met, D2 met (113 of 134 confirmable reads, 84%), the two small clusters stay NOVEL** (no alignment on the primary or on either haplotype assembly; supported consensus, read divergence 0.0008 and 0.0009; 2,445 and 2,330 bp). **Clean control seed 3: no CONFIRMED, NOVEL or UNSUPPORTED (bar met), 97% of reads PRESENT or ALLELE-LIKE. Fresh clean control seed 4: bar NOT met** (one CONFIRMED cluster of 9 reads, 1.3% of the 684 clustered reads; 98.7% of the reads are PRESENT or ALLELE-LIKE, that clause met). The false CONFIRMED is the diverged ortholog again: chr11 (mother 29.19 Mb, father 23.12 Mb), 1,977 bp, the consensus is the father's haplotype (Tp 1.0000) and the primary carries the mother's version at 0.9889, 1.1% away, just beyond the allele cutoff, so "absent from the primary, present in another assembly of the same animal" is true and it is an allele. **Verdict: the discovery run finds the real reference-absent sequence and two sequences no assembly holds, and its false-discovery rate on reads that are clean on the reference is about 1% of the clustered reads, all of it at the 1% allele boundary; the registered specificity bar is not met on the fresh control.**

---

### Amendment 32 (2026-10-08, before any k-mer count): DNA verification of the discovery candidates

**Question.** Are the sequences reconstructed from unmapped reads present in the individual's own DNA? Data: the whole-genome reads of KB3781 ("Jim", SRR26039725, paired Illumina, 138.8 M spots, 41.9 Gb, 151 bp; k-mer depth per haploid genome lambda = 4.98, Amendment-free measurement of 2026-09-25 on 516 single-copy autosomal families). Method: the canonical 21-mers of each sequence are looked up in a streamed pass over all the reads (`kmerhit.c`, exact match, both strands, count per k-mer); no alignment.
**Sequences.** T1 = the CONFIRMED 113-read consensus (5,694 bp; positive control with independent evidence: identical to the mother's chr12 locus); N1 and N2 = the two NOVEL consensus sequences (2,445 and 2,330 bp); P1 = a PRESENT consensus of the clean control (the 61-read, 2,393 bp cluster; positive control for a sequence the primary holds); S1 and S2 = dinucleotide-shuffled copies of T1 and N1 (same composition, negative controls).
**Statistic (robust to splicing).** A cDNA consensus contains k-mers that span exon junctions and exist in no genomic read. A base is COVERED if at least one of the 21-mers that contain it has count >= 2 (count 1 is read-error range); the covered fraction of the sequence and the median count of the k-mers with count >= 2 are reported.
**Bars (set now).** (V1) T1: covered >= 95% and median count within [2.5, 25] (0.5 to 5 times lambda: one haplotype to a multi-copy locus). (V2) the shuffled controls: covered <= 2%. (V3) P1: covered >= 95%. **Classification of the candidates (no bar, a reading):** DNA-SUPPORTED = covered >= 90% and median count in [2.5, 25]; PARTIAL = covered in [10%, 90%); DNA-ABSENT = covered < 10%.
**Decision.** The verification is reported as working iff V1 to V3 hold; then N1 and N2 are classified. A DNA-ABSENT candidate is an RNA-only sequence (contamination, a trans-spliced or chimeric read product); DNA-SUPPORTED is a sequence of the animal's genome absent from all three assemblies.

#### Amendment 32 result (2026-10-08; `kmerhit.c`, `dna_verify.py`; logs `dna_verify.jim.log`, `candidates_vs_human.log`)

All 277,646,976 reads (41.9 Gb) of KB3781's whole-genome run were streamed once (6 chunks, 2 x 2 min); the 20,674 distinct canonical 21-mers of the six sequences were counted exactly.

| sequence | length | k-mers with count >= 2 | covered bases (count >= 2) | median count of present k-mers |
|---|---|---|---|---|
| **T1** (113-read CONFIRMED) | 5,694 | **97.1%** | **99.8%** | **6.0** (one haplotype; lambda 4.98) |
| **P1** (PRESENT control) | 2,393 | 96.0% | 99.9% | 11 (two haplotypes) |
| **N1** (NOVEL, 4 reads) | 2,445 | **0.3%** | **5.2%** | 8 (7 k-mers) |
| **N2** (NOVEL, 3 reads) | 2,330 | **0.2%** | **3.7%** | 11 (5 k-mers) |
| S1 (shuffled T1) | 5,694 | 0.3% | 4.8% | 9 |
| S2 (shuffled N1) | 2,445 | 0.1% | 1.7% | 8.5 |

**V1 met** (T1: 99.8% covered, median 6.0 in [2.5, 25]). **V3 met** (P1 99.9%). **V2 not met as written**: the shuffled negatives are covered 4.8% (S1) and 1.7% (S2) against the bar <= 2%: a handful of random 21-mers (0.1 to 0.3% of them, chance matches in a 3.5 Gb genome read at 12x) each cover 21 bases, which the covered-bases statistic magnifies; on the k-mer fraction the separation is 97% against 0.3% and the bar's intent holds. **Reading: N1 and N2 are DNA-ABSENT**: their covered fractions (5.2%, 3.7%) and k-mer fractions (0.3%, 0.2%) equal the shuffled negatives'. T1 aligns to human chr2 (98.2%, an ortholog); N1 and N2 do not align to the human genome either (no record), have ordinary composition (GC 0.53 and 0.58, no homopolymer or low-complexity run, no polyA), so Illumina GC dropout and human contamination are not the explanation.
**What the discovery run therefore found, graded by independent evidence.** One sequence (113 reads, 5.7 kb) absent from the primary assembly, identical to the mother's assembly and present in the animal's DNA at one-haplotype depth: a real reference-absent sequence reached from the transcriptome. Two small clusters (7 reads) that are in no assembly, no human sequence and none of the animal's 41.9 Gb of DNA: RNA-only products (contamination of the culture or the library, or chimeric or trans-spliced read products), not unannotated gorilla genes; the verification step removes them as candidates. The verification itself cost 4 minutes of streaming and a 17 KB k-mer set.

---

### Amendment 33 (2026-10-08, before any result): the discovery run on the other animals

**Data.** Iso-Seq reads of other animals, aligned by the project with `splice:hq -uf` to their own reference, with the primary-unmapped records kept: gorilla (the testis animal, not KB3781; `GGO_mm.bam` against the same mGorGor1 primary) 5,506 unmapped reads; human (A119b against CHM13, `A119b.t2t.bam`) 9,841; orangutan (`PPY_mm.bam`, mPonPyg2) 23,051; chimpanzee (`PTR_mm.bam`, mPanTro3) 449. (KB3781 itself was Amendment 28.)
**Pipeline.** Exactly the adopted one: the unmapped reads, `--proper` edge clustering with the star step, abPOA consensus, the support requirement, terminal-exon rescue, the 5' G rule with the library gate computed per animal from 200,000 of its mapped primaries. Each consensus is aligned to the animal's own reference (R) and to every other assembly on disk (the same species' other individual, and the other primates): E = the best identity x coverage over the others. Classes: PRESENT (R >= 0.999), ALLELE-LIKE (R in [0.99042, 0.999)), ELSEWHERE (R < 0.99042 and E >= 0.9: absent from the animal's own assembly, a homolog or the same sequence in another assembly), NOVEL (R and E both < 0.9 or no hit), UNSUPPORTED (Amendment 21).
**Controls (specificity), per animal.** 959 reads sampled (seed 5) from the animal's mapped primaries with divergence <= 0.00958, same pipeline: at least 95% of the clustered reads must be PRESENT or ALLELE-LIKE.
**Prediction (one, from Amendment 32).** The STON1-GTF2A1L readthrough that is absent from the mGorGor1 primary (5.7 kb, 113 reads in KB3781) is absent for every animal that shares that reference: the testis gorilla's unmapped reads contain it (a supported cluster of >= 20 reads that is ELSEWHERE against KB3781's mother's assembly at >= 0.99).
**Bars.** (E1) the prediction above. (E2) each animal's control: >= 95% of clustered reads PRESENT or ALLELE-LIKE. Reported, not judged: for each animal the number of clusters and reads in each class, the largest ELSEWHERE and NOVEL clusters with their best hits (assembly, locus), the share of unmapped reads that end in a cluster; for the gorilla testis animal, the DNA check of every NOVEL and ELSEWHERE cluster against KB3781's whole-genome reads (a different individual: presence supports, absence is not conclusive).
**Decision.** The discovery path is reported as transferring iff E1 and E2 hold for the animals run; nothing is called novel without the DNA or cross-species evidence stated.

#### Amendment 33b (2026-10-08, before any result): implementation choices fixed before the run

- **Reads.** Unmapped primaries (`-f 4 -F 2304`) of `GGO_mm.bam` (5,506; median 2.9 kb), `A119b.t2t.bam` (9,841; median 72 bp), `PTR_mm.bam` (449; median 2.1 kb), `PPY_mm.bam` (23,051; median 68 bp). Nothing is length-filtered: short reads simply do not form clusters of >= 3.
- **Gate (per animal).** The 5' signature (`libsig`) of the first 200,000 clean mapped primaries (MAPQ >= 10, flag & 2308 = 0, de <= 0.0096) in 5 Mb windows centred on each of the 12 longest contigs. `discover.py prepare <animal>`.
- **Control (per animal).** From the same windows, clean reads (de <= 0.00958) lying in a 5 kb bin with >= 5 such reads (so they can cluster), 959 sampled with seed 5 (all if fewer). Same pipeline and classes. It is clean on the animal's own reference only (no haplotype assemblies exist for PTR, PPY, A119b). The easier-than-real caveat of the earlier controls applies.
- **Others (E).** gorilla testis animal: Tm and Tp, the maternal and paternal haplotype assemblies used in Amendment 28 (the "mother's assembly" of the E1 prediction is Tm; the testis animal is not KB3781). Human A119b: the mGorGor1 primary. Chimpanzee and orangutan: CHM13 and the mGorGor1 primary. Orangutan and chimpanzee reference indexes are built from the FASTA at alignment time (no 13 GB index kept; disk is the limit). Cross-species scores are homology, not identity of the same locus: ELSEWHERE against a different species means "a homolog exists", and that reading is the one reported.
- **Classes.** `flagmetric.elsewhere_class`: PRESENT (R >= 0.999), ALLELE-LIKE (R in [0.99042, 0.999)), ELSEWHERE (E >= 0.9), NOVEL (R and E < 0.9), DIVERGED (R in [0.9, 0.99042) and E < 0.9: a distant paralog or an incomplete alignment inside its own reference; the class existed before and Amendment 33 left it unnamed), UNSUPPORTED (median read divergence > 0.00958 or none). Unit tests in `test_flagmetric.py`, `test_discover.py`.
- **E1 (made exact).** In `ggo_testis`: a supported cluster of >= 20 reads with Tm >= 0.99 and R < 0.99042.

#### Amendment 33 result (2026-10-08; `discover.py`, logs `docs/unmapped_rescue_results/discover_animals.*.log`)

| run | reads -> clustered | clusters | ELSEWHERE (cl/reads) | NOVEL | UNSUPPORTED | PRESENT + ALLELE-LIKE share |
|---|---|---|---|---|---|---|
| gorilla testis | 5,506 -> 61 (1.1%) | 9 | 6 / 51 | 3 / 10 | 0 | n/a |
| gorilla testis control | 959 -> 659 | 79 | 3 / 22 | 0 | 0 | **96.7% (E2 met)** |
| human A119b | 9,841 -> 742 (7.5%) | 78 | 56 / 584 | 7 / 35 | 15 / 123 | n/a |
| human control | 959 -> 322 | 55 | 5 / 18 | 1 / 5 | 5 / 20 | **86.7% (E2 NOT met)** |
| chimpanzee | 449 -> 58 (12.9%) | 7 | 5 / 26 | 1 / 8 | 1 / 24 | n/a |
| chimpanzee control | 959 -> 495 | 78 | 1 / 5 | 2 / 30 | 4 / 21 | **88.7% (E2 NOT met)** |
| orangutan | 23,051 -> 53 (0.2%) | 12 | 7 / 32 | 2 / 8 | 3 / 13 | n/a |
| orangutan control | 959 -> 324 | 49 | 2 / 19 | 0 | 0 | **94.1% (E2 NOT met, 0.9 point short)** |

**E1 met.** The testis gorilla's largest ELSEWHERE cluster has 30 reads, 5,523 bp, no alignment on the primary (R none), Tm = 0.998 at chr12_mat 95,963,504-96,028,675: the STON1-GTF2A1L readthrough locus of Amendment 32. **E2 met in 1 of 4 animals.** The failing controls are not "specific": in every control the non-PRESENT clusters are consensus of 3 to 6 reads from repeat-rich multi-copy loci (chr6, chr7, chrX in the human control) whose consensus aligns to the SAME locus on its own reference at 0.68 to 0.99 (the support rule catches the worst); the rule `ELSEWHERE iff E >= 0.9` fires on them when a second assembly is only slightly worse than R (R 0.988, E 0.965). This is a defect of the rule as written (ELSEWHERE should mean "better explained elsewhere than on R"), disclosed, not repaired here.
**Observation (after the numbers, not a bar).** Across the 4 controls (261 clusters) 0 ELSEWHERE or NOVEL clusters have NO alignment on the own reference; across the 4 real runs 86 of 87 do.
**The unmapped sets are a few hundred reads each.** Most unmapped records are short (human 72 bp, orangutan 68 bp median); 0.2% to 13% of them form clusters.
**Human caveat found while reading.** The largest human clusters (146 reads, DNAJC15 281 bp with 35 reads) are genes every human has, at 0.96 to 0.99 of the gorilla primary and no alignment on CHM13 with `splice:hq`. With the lenient preset `splice`, 19 of the 78 human consensus sequences align to CHM13 (4 of 78 with `splice:hq`). So "absent from the reference" in these libraries is partly the aligner setting (`splice:hq` needs ~99% identity) and partly sequence that is not human.

### Amendment 34 (2026-10-08, before any number): whose reads are the unmapped clusters?

**Question.** Are the ELSEWHERE / NOVEL / UNSUPPORTED clusters of Amendment 33 sequence of the animal's own species (absent from its reference, which is what a discovery would be), or reads of another primate (a library property, not a discovery)? A human read hits CHM13 at >= 0.999 under any preset, so a human cluster with no CHM13 alignment is not a human sequence.
**Test.** The consensus of every ELSEWHERE, NOVEL and UNSUPPORTED cluster of the four real runs (pooled, prefixed by the animal), aligned with the lenient preset `minimap2 -c -x splice -uf -N 5` to CHM13 (HSA), the mGorGor1 primary (GGO), the maternal haplotype assembly (Tm, counted as GGO), the chimpanzee mPanTro3 (PTR) and the orangutan mPonPyg2 (PPY). Score of a cluster on a species = (sum of matching bases over the primary records, `tp:A:P`) / consensus length, capped at 1 (approximate identity x coverage, no end rescue); a species' score is the best over its assemblies. Best species = the highest score. **Same-species** = best species is the animal's own and its score >= 0.99.
**Predictions (made now).** (P1) gorilla testis: >= 80% of the ELSEWHERE reads lie in same-species clusters (they are the maternal-haplotype sequence: Tm >= 0.998 already). (P2) human: >= 50% of the ELSEWHERE + NOVEL reads lie in clusters whose best species is NOT human with score >= 0.99 (the reads are another primate's). Reported without a prediction: chimpanzee and orangutan; the species composition of the NOVEL clusters; the clusters whose best score is under 0.9 on every species.
**Decision.** Where P2 holds, the human ELSEWHERE clusters are not discoveries and the discovery path is reported with that contamination class named; where P1 fails, the gorilla ELSEWHERE count is re-read.

#### Amendment 34 result (2026-10-08; `species_check.py`; `docs/unmapped_rescue_results/species_check.{log,json}`)

106 consensus sequences (ELSEWHERE, NOVEL, UNSUPPORTED of the four real runs) aligned with `splice -uf` to CHM13, the mGorGor1 primary and maternal haplotype (GGO), mPanTro3, mPonPyg2.
**P1 met:** 47 of 51 gorilla-testis ELSEWHERE reads (92%) lie in clusters whose best species is gorilla at >= 0.99 (the other 4 reads: 0.9 to 0.99). The 30-read cluster aligns end to end to T1, the DNA-confirmed 113-read consensus of KB3781 (Amendment 32): 5,523 of 5,523 bases at **99.82%** identity. The same reference-absent readthrough is expressed in two gorillas (different animals, different tissues), neither of which the primary carries.
**P2 NOT met (0.155, bar 0.50), and the reading is different from the one predicted.** The human ELSEWHERE + NOVEL reads (619) are not a single other species at >= 0.99: 507 reads (82%) fall in the 0.9 to 0.99 band on their best species, 85 reads (14%) are chimpanzee >= 0.99 (the 35-read 281 bp DNAJC15 cluster is 0.996 on chimpanzee, a short sequence), 11 gorilla, 16 under 0.9 everywhere. **CHM13 scores 0.0 on all the top clusters** (146 reads: HSA 0.0, GGO 0.986, PPY 0.943, PTR 0.777 partial). Sequence of genes every human carries that no human assembly hits and that sits at 0.96 to 0.99 of three apes, with no panel species at the 0.999 a same-species read would give, is not a human discovery; the panel cannot name its species (a primate outside it, or a mixed-species library). **Unresolved.** Next step that would settle it: BLAST of the 3 kb top consensus against nt.
**The other animals' clusters:** chimpanzee, orangutan: no cluster is same-species (0 of 6 and 0 of 9 ELSEWHERE + NOVEL; their own reference holds none of them, and their best species is another ape at 0.9 to 0.99). The gorilla testis NOVEL clusters (3, 10 reads: 269, 280, 1,499 bp) match no gorilla assembly and sit at 0.9 to 0.97 on chimpanzee / orangutan: the same class.
**What survives for the thesis.** The reference-absent sequence found in the 959 KB3781 reads (Amendment 32) is reproduced in an independent gorilla at 99.82% over its full length, and 6 of 6 of that animal's ELSEWHERE clusters are the maternal haplotype's sequence, as predicted. In the other three species the discovery path returns no same-species sequence that the animal's own reference lacks: what it finds in unmapped reads is mostly other species' sequence. E2 (specificity of the control) was met only in the gorilla; the rule `ELSEWHERE iff E >= 0.9` should also require E > R (defect found in the controls, not repaired).

### Amendment 35 (2026-10-08, before any re-scoring): ELSEWHERE means "explained better elsewhere than on R"

**Defect (Amendment 33 result).** `ELSEWHERE iff R < 1 - delta and E >= 0.9` fired on clean-control consensus sequences whose own reference holds them nearly as well as another assembly (R 0.988, E 0.965; R 0.974, 0.982, 0.962 with E within 0.001 of R): 11 clusters in the four controls (3, 5, 1, 2). These were seen before this rule was written, so F1 below is a bug-fix check, not a blind test; F3 is the blind one.
**Rule.** ELSEWHERE iff R < 1 - delta and E >= 0.9 **and E > R** (a missing score is 0; no new threshold). Then DIVERGED iff R >= 0.9 and not ELSEWHERE (the own reference holds a hit that no other assembly beats: paralog, structural difference, consensus of mixed copies); NOVEL iff R < 0.9 and E < 0.9. PRESENT, ALLELE-LIKE, UNSUPPORTED unchanged.
**Bars (set now).** (F1) re-scoring the saved scores of the four controls (seed 5): 0 ELSEWHERE clusters in each. (F2) re-scoring the four real runs: only clusters with R >= 0.9 change class; the gorilla-testis E1 cluster (30 reads, R none) stays ELSEWHERE; the list of changed clusters is reported. (F3, blind) fresh controls, seed 6, one per animal, same pipeline: 0 ELSEWHERE clusters in each. **Not expected to change and not claimed:** the E2 share (PRESENT + ALLELE-LIKE >= 95%), which this fix cannot move; it is reported again for the seed-6 controls.
**Decision.** The rule replaces Amendment 33's iff F1, F2 and F3 hold; otherwise it is reported as not adopted.

#### Amendment 35 F1 result and Amendment 35b (2026-10-08, after F1 on the saved scores, before F2/F3 and before any seed-6 control exists)

**F1 NOT met by the strict rule E > R.** Re-scoring the saved controls: ELSEWHERE clusters 1 (gorilla testis), 0 (human), 0 (chimp), 0 (orangutan); the 11 clusters became DIVERGED except one, 8 reads / 1,424 bp in the gorilla control, R = 0.982, Tp = 0.983: E exceeds R by 0.001, a near-tie that sampling noise decides.
**Rule 35b (post hoc to that F1 result, so the blind F3 is what counts).** ELSEWHERE iff R < 1 - delta and E >= 0.9 and (E >= 0.999 or E > R + delta): another assembly holds it completely, or beats the own reference by more than the registered allele tolerance (0.00958), the same tolerance that separates ALLELE-LIKE from PRESENT; two registered constants, no new one. DIVERGED / NOVEL as in Amendment 35.
**Bars unchanged and applied to 35b:** F1 (saved controls, 0 ELSEWHERE each; not blind), F2 (real runs: only clusters with R >= 0.9 change, E1 cluster stays ELSEWHERE), F3 (fresh seed-6 controls, one per animal, 0 ELSEWHERE each; blind). Decision as in Amendment 35.

#### Amendment 35 F2/F3 result for 35b and Amendment 35c (2026-10-08, after the seed-6 numbers, before any seed-7 control exists)

**35b results.** F1 met (0 ELSEWHERE in the four saved seed-5 controls). F2 met (no real-run cluster changes class; the one human cluster with an own-reference hit has R = 0.062, E = 0.989 and stays ELSEWHERE; E1 cluster unchanged). **F3 NOT met: fresh seed-6 controls give ELSEWHERE 0 (gorilla, 637 clustered reads), 0 (human, 296), 2 clusters / 9 reads (chimp, 521), 0 (orangutan, 325).** The two chimp clusters hit the chimp reference at 0.971 and 0.970 and a HUMAN / GORILLA assembly higher (0.968 / 0.981 and 0.994 / 0.994): a cross-species homolog can outscore the animal's own reference through coverage, so for another species E > R + delta is not evidence of absence. Under the Amendment 35 decision rule, 35b is **not adopted**. (Share PRESENT + ALLELE-LIKE of the seed-6 controls: gorilla 97.96%, human 93.6%, chimp 87.7%, orangutan 92.3%.)
**Rule 35c (post hoc to the seed-6 result; the seed-7 controls are the blind test).** The other assemblies split by species. *Same species* (only the gorilla testis animal has them: Tm, Tp): the 35b rule, unchanged (E_same >= 0.9 and (E_same >= 0.999 or E_same > R + delta)); it keeps the case the thesis cares about, a copy absent from the primary with a divergent paralog on it and the whole sequence on the other haplotype. *Another species* (human, gorilla primary for chimp / orangutan / human; configured per animal): ELSEWHERE only if the own reference does not hold the sequence at all, R < 0.9 (the constant that already separates NOVEL from DIVERGED), and E_cross >= 0.9. DIVERGED iff R >= 0.9 and not ELSEWHERE; NOVEL iff R < 0.9 and neither holds >= 0.9.
**Bars for 35c (set now).** (F1') the saved controls (seeds 5 and 6; not blind): 0 ELSEWHERE in all eight. (F2') the real runs: no cluster changes class relative to Amendment 33 except those with R >= 0.9 (none expected: human has one with R = 0.062). (F3', blind) fresh controls, seed 7, one per animal: 0 ELSEWHERE clusters in each. Decision as in Amendment 35: adopted iff F1', F2', F3' all hold.

#### Amendment 35c result (2026-10-08; `flagmetric.elsewhere_class`, `discover.py rescore`; `docs/unmapped_rescue_results/discover_animals.rescore_35c.log`)

**F1' met.** The eight saved controls (seeds 5 and 6): 0 ELSEWHERE clusters in each (was 3, 5, 1, 2 and 0, 0, 2, 0 under the Amendment 33 rule and, for seed 6, under 35b the chimp's 2). **F2' met.** The four real runs: no cluster changes class (gorilla testis 6 ELSEWHERE / 51 reads, human 56 / 584, chimp 5 / 26, orangutan 7 / 32 as before; the E1 cluster stays ELSEWHERE). **F3' (blind, fresh seed-7 controls) met: ELSEWHERE 0 of 87 clusters (gorilla, 638 clustered reads), 0 of 51 (human, 287), 0 of 85 (chimp, 530), 0 of 62 (orangutan, 355).** **Rule 35c is adopted** (it is what `discover.py classify` now computes; the Amendment 33 table is the superseded rule).
**What the fix did not do, and one thing it exposed.** (i) The E2 share (PRESENT + ALLELE-LIKE >= 95%) is unchanged by construction. Over seeds 5 / 6 / 7: gorilla 96.7 / 98.0 / 98.3% (met in all three); human 86.7 / 93.6 / 92.3%; chimp 88.7 / 87.7 / 89.1%; orangutan 94.1 / 92.3 / 85.4% (never met). The shortfall is DIVERGED, UNSUPPORTED and NOVEL clusters of clean reads, i.e. consensus of a few reads from multi-copy loci. (ii) **NOVEL is not specific either:** clean controls contain NOVEL clusters in human (1, 1, 1 cluster; 5, 4, 3 reads), chimp (2, 1, 3; 30, 3, 15), orangutan (0, 0, 2; 0, 0, 11), gorilla (0, 0, 0). So the NOVEL clusters of the real runs (human 7, chimp 1, orangutan 2, gorilla testis 3) are not evidence of novel sequence on their own; the gorilla class is the only one with a clean control record.

### Amendment 36 (2026-10-08, before any number): why do clean-read control clusters have a bad consensus?

**Problem (Amendments 33, 35c).** In the 12 control runs (4 animals x seeds 5, 6, 7; reads clean on their own reference, each mapped uniquely at MAPQ >= 10) the share of clustered reads in PRESENT or ALLELE-LIKE clusters is 85 to 98%; the rest lie in DIVERGED, UNSUPPORTED and NOVEL clusters, which is what makes NOVEL non-specific.
**Hypotheses.** H1 (mixed copies): such a cluster joins reads whose own alignments lie at different loci (paralogs within the edge tolerance delta = 0.00958 of each other), so its consensus is a blend that matches no single copy. H2 (consensus quality): the reads of the cluster map to ONE locus and the consensus is still poor (abPOA, ends, coverage).
**Diagnosis (descriptive, no repair yet).** For each control cluster, the primary alignment of every read (contig, position) is recovered from the BAM. A cluster is *multi-locus* if its reads' alignments are not all within 100 kb on one contig. **D1:** among reads in DIVERGED + UNSUPPORTED + NOVEL clusters, >= 50% lie in multi-locus clusters, versus < 10% among reads in PRESENT + ALLELE-LIKE clusters. If D1 holds, H1 is the main cause; if the failing clusters are mostly single-locus, H2. Reported: per class the number of clusters, reads and multi-locus share; the distribution of the number of distinct loci; for single-locus failing clusters the consensus score and median read divergence.
**Next step (registered after the diagnosis, not now).** If H1: split a multi-locus cluster by its reads' alignments' differing columns (the PART machinery of Amendments 18 to 20, guarded by block purity) and re-classify the leaves; bars set after D1.

#### Amendment 36 result (2026-10-08; `control_purity.py`, `test_purity.py`; probes on the 12 control runs; anatomy in `/mnt/linuxdisk/tmp/o3_rescue/{control_purity,control_anatomy,unsupported_anatomy}.json`)

**D1 NOT met; H1 (mixed copies) refuted.** Multi-locus read share: DIVERGED + UNSUPPORTED + NOVEL **8.9%** (2 of 72 clusters, 404 reads), PRESENT + ALLELE-LIKE 2.9% (27 of 745 clusters, 4,985 reads). 70 of the 72 failing clusters are single-locus. The failures are not mixtures of paralogous loci.
**What they are (descriptive).** (a) Not coverage: the failing consensus sequences align over 99.5%+ of their length; their aligned identity is the loss (DIVERGED 0.979, UNSUPPORTED 0.83, NOVEL 0.88; PRESENT 0.9997). (b) **DIVERGED (32 clusters, 173 reads):** the consensus is supported by its reads (median read divergence 0.0018) and 22 of 32 have gap-compressed divergence `de` <= 0.00958 on R, the unit the allele cutoff is defined in, while the per-base score (mismatches + every gap base, 0.0206 median) calls them diverged: a few indel events of tens of bases cost 1 event in `de` and 1 to 2% in the score. A units mismatch between the cutoff and the score, not a consensus failure. (c) **NOVEL (11 clusters, 71 reads):** 9 of 11 are supported by their reads, aligned identity 0.88: large indels relative to R; same units issue, with insertions. (d) **UNSUPPORTED (29 clusters, 160 reads): the consensus itself.** The abPOA consensus lies within delta of 28 of the 160 reads (17%); for every one of the 29 clusters a single read of the cluster lies within delta of every other read (100%; true by the star construction, so this bounds what a consensus can reach). The reads agree with each other; the consensus does not agree with them (12 reads all at 0.063 from it).
### Amendment 37 (2026-10-08, before any consensus is rebuilt): the consensus mode

**Hypothesis H3.** `abpoa_consensus` aligns every read globally (`aln_mode="g"`) to a graph seeded by the longest read. Iso-Seq reads of a locus are fragments of unequal extent (read lengths in one cluster 924 to 3,368 bp), and a global alignment forces end gaps on the fragments, which corrupts the heaviest-bundle consensus. Local (`l`) or extension (`e`) mode should not.
**Test, in order.** (S1) a unit test on a synthetic fragment cluster (truth transcript, reads of unequal extent, 0.5% errors, fixed seeds): global-mode consensus identity to the truth < 0.99, the best other mode >= 0.995. If global is already >= 0.99 there, H3 is not reproduced and the real-data test is not run. (S2) the 29 UNSUPPORTED clusters: reads within delta of the consensus, global 28 of 160 (17%) -> **>= 80% with the chosen mode**, and >= 24 of 29 clusters supported (median read divergence <= 0.00958). (S3) no harm: 200 PRESENT / ALLELE-LIKE control clusters (seed 1): reads within delta of the consensus not lower than with global mode minus 0.01, and their identity x coverage on R (re-aligned) median not lower by more than 0.0005. (S4) the 12 control runs re-classified with the chosen mode: the PRESENT + ALLELE-LIKE share does not fall in any run; the 30-read testis cluster stays ELSEWHERE (Tm >= 0.99); the real-run class counts reported. (S5, blind) fresh seed-8 controls, one per animal, both modes on the same reads: new share >= old in 4 of 4.
**Chosen mode** = the mode with the largest S2 reads-within-delta among `l` and `e` (ties: `e`), fixed before S3 to S5. **Decision:** adopted as the default consensus iff S1 to S5 hold. **Not addressed here (reported):** the units mismatch in DIVERGED and NOVEL (b, c); it changes the meaning of PRESENT / ALLELE-LIKE and needs its own amendment with a sensitivity check on erased copies.

#### Amendment 37 S1 result and 37b (2026-10-08, after S1 on 8 seeds, before the 200-seed run)

**S1 as written NOT met.** 8 synthetic fragment clusters (seeds 1 to 8; 12 reads, 2,000 bp, starts 0 to 500 bp late, ends up to 200 bp early, 0.5% errors): median identity to the truth g 0.9989, l 1.0000, e 1.0000. By the Amendment 37 rule the real-data test is not run on this S1. But the worst case differs: **g 0.7603 (one of 8 clusters garbled)**, l 0.9995, e 0.9918. The real clusters are the same shape: 29 of 817 (3.5%) garbled, the rest fine. The defect is a tail, and a median cannot see it; S1 was the wrong statistic.
**37b (replaces S1; the 8-seed numbers are not blind, the new run is).** S1': 200 fresh synthetic fragment clusters (seeds 101 to 300, same generator, same parameters). Failure = identity to the truth < 0.99. Bars: global failure rate >= 3% (the rate the real clusters show) and **the better of l / e fails in <= 1%**. If S1' holds, S2 to S5 are run as registered; if not, the consensus-mode hypothesis is reported as not supported.

#### Amendment 37 result (2026-10-08; `partition.abpoa_consensus(mode=)`, `consensus_modes_probe.py`, `reclassify_batch.py`, `discover.py`; logs `docs/unmapped_rescue_results/discover_animals.reclass_local*.log`)

**All bars met; abPOA local mode is adopted for the consensus of `discover.py` (`CONS_MODE = "l"`).** The frozen global mode stays the default of `partition.py`, `consensus.py` and the bed / synthetic pipelines (their reads are full length; not re-validated, not changed).
- **S1' (200 fresh synthetic fragment clusters): met.** Failure (< 0.99 identity to the truth) g 26% (52 of 200), l 0% (0), e 13% (26). The earlier 8-seed median hid the tail (S1 as first written was not met and was replaced by 37b before this run).
- **S2 (29 UNSUPPORTED control clusters): met.** Reads within delta of the consensus: g 28 of 160 (17%), **l 156 (97.5%)**, e 144 (90%); clusters supported (median <= delta): g 0, **l 29 of 29**, e 27. Mode chosen by the registered rule: l.
- **S3 (no harm, 200 good clusters + all 745): met.** Reads within delta g 0.998, l 0.999. Identity x coverage on R of the 745 clusters that are good under both consensus sequences: median change +0.00000. **A regression, reported:** 3 of 745 good clusters get a worse score (ptr_control_s7: 1,244 -> 2,650 bp, ALLELE-LIKE -> NOVEL; ppy_control: 1,704 -> 2,170 bp, ALLELE-LIKE -> DIVERGED; one 0.9984 -> 0.9922): local mode keeps the unsupported tail of the longest read (870 of its 2,650 bases do not align to R; the "clean" control only requires the aligned part to be clean). The mean change is -0.00018 (worst -0.14).
- **S4 (12 control runs re-classified): met.** PRESENT + ALLELE-LIKE share, global -> local, seeds 5 / 6 / 7: gorilla 0.967 / 0.980 / 0.983 -> **1.000 / 1.000 / 1.000**; human 0.866 / 0.936 / 0.923 -> **1.000 / 1.000 / 1.000**; chimp 0.887 / 0.877 / 0.891 -> **0.990 / 0.973 / 0.968**; orangutan 0.941 / 0.923 / 0.854 -> **0.985 / 0.991 / 1.000**; no run falls. The testis E1 cluster stays ELSEWHERE (30 reads, 5,520 bp, R none, Tm 0.998). 67 clusters become good, 2 stop being good.
- **S5 (blind, fresh seed-8 controls, same reads, both modes): met 4 of 4.** gorilla 0.995 -> 1.000; human 0.926 -> 0.991; chimp 0.899 -> 0.971; orangutan 0.933 -> 0.992. **E2 (>= 95%) is now met in all 16 control runs** (it was met in 1 of 4 animals).
**Correction to Amendment 36 (b), (c).** DIVERGED and NOVEL clusters of the controls were also garbled global-mode consensus, not a units mismatch: with the local consensus they disappear (NOVEL in controls: 2 clusters / 12 reads in 16 runs, was 11 clusters / 71 reads in 12; DIVERGED 0 to 3 clusters per run). The `de` versus per-base remark stays true for the cases that remain but is not needed to reach the bar.
**Real runs with the local consensus (class counts, clusters / reads):** gorilla testis ELSEWHERE 6 / 51, NOVEL 3 / 10 (unchanged); human ELSEWHERE 56 / 584 -> **70 / 703**, NOVEL 7 / 36, DIVERGED 1 / 3, UNSUPPORTED 15 / 123 -> **0**; chimp ELSEWHERE 5 / 26 -> 6 / 50, NOVEL 1 / 8, UNSUPPORTED 1 / 24 -> 0; orangutan ELSEWHERE 7 / 32 -> 11 / 49, NOVEL 2 / 8 -> 1 / 4, UNSUPPORTED 3 / 13 -> 0. The UNSUPPORTED clusters were real clusters whose consensus was wrong; they now classify as ELSEWHERE (no hit on the own reference). **The Amendment 34 species check used the global consensus and is not yet refreshed**; its human reading (no CHM13 hit, no panel species at >= 0.99) is expected to stand and is unverified.
**Open.** (i) End support: trim the consensus to columns covered by >= 2 reads (the 3 regressions). (ii) Refresh Amendment 34 on the local consensus. (iii) The same global-mode defect probably exists in the bed / famsim pipelines when reads are fragments (not tested).

### Amendment 38 (2026-10-08, before any trimming code runs): trim the consensus to columns covered by at least two reads

**Defect (Amendment 37 S3).** abPOA local mode keeps the unsupported flank of the seed (longest) read: 3 of 745 good control clusters got a longer consensus (1,244 -> 2,650 bp; 1,704 -> 2,170 bp; 1,833 -> 4,345 bp) whose extra bases align nowhere on R and drop the class (ALLELE-LIKE -> NOVEL / DIVERGED); in two of them the longest read has a 450 to 870 bp tail that does not align.
**Rule.** After the local-mode consensus, align the cluster's reads (primary alignments, `splice:hq -uf`, the same mapping as the support test) and trim both ends of the consensus to the first and last column covered by >= 2 reads (reference span of each alignment, M / = / X / D / N). The window is empty or the cluster has < 2 aligned reads: keep the consensus. One replicate is the minimum that separates a molecule from a sequence; no other threshold. Variant name `lt` (local + trim).
**Bars (set now).** (T1, synthetic, 200 fresh seeds 301 to 500, the Amendment 37b generator plus one read that carries a 300 bp random unsupported tail on a random side): identity to the truth < 0.99 counts as failure; **untrimmed local mode fails in >= 50%, the trimmed consensus in <= 2%.** (T2, no harm, the 200 tail-free seeds 101 to 300): trimmed failure rate <= 1%; the median consensus length loss relative to untrimmed is reported. (T3, the 16 control runs re-classified, seeds 5 to 8): (a) the 3 regressions return to PRESENT / ALLELE-LIKE; (b) the PRESENT + ALLELE-LIKE share is not lower than with local mode in any run; (c) among clusters good under local mode, **none** stops being good, and the median change of identity x coverage on R is >= -0.0005. (T4, blind, fresh seed-9 controls, one per animal, both consensus variants on the same reads): trimmed share >= untrimmed in 4 of 4. Reported, no bar: the real-run class counts; the E1 cluster must stay ELSEWHERE with Tm >= 0.99; the number of real ELSEWHERE / NOVEL clusters shortened by more than 5%; the number of clusters whose trim is empty.
**Decision.** Trimming is added to `discover.py` iff T1 to T4 hold.

#### Amendment 38 result (2026-10-08; `endtrim.py`, `test_endtrim.py`, variant `lt` of `reclassify_batch.py`; logs `discover_animals.reclass_lt_*.log`): NOT adopted

- **T1 met.** 200 fresh synthetic fragment clusters (seeds 301 to 500, the longest read carries a 300 bp random tail on a random side): identity to the truth < 0.99 in **63.5%** (127) of untrimmed local-mode consensus sequences, **0%** of trimmed; median identity 0.894 -> 1.000.
- **T2 met.** 200 tail-free seeds (101 to 300): trimmed failure 0%; median consensus length loss **1.8%** (mean 2.5%, max 9.7%).
- **T3 NOT met.** (a) Of the two class-changing local-mode regressions, one recovers (ppy_control 2,170 -> 1,234 bp, DIVERGED 0.931 -> ALLELE-LIKE 0.998); the other does not (ptr_control_s7 2,650 -> 2,563 bp, still NOVEL, R 0.853: at least two reads cover the extra sequence). (b) The PRESENT + ALLELE-LIKE share falls in **3 of 12 runs** (human 1.000 -> 0.991 and 1.000 -> 0.983; orangutan seed 6 0.991 -> 0.978), is equal in 7, rises in 1 (orangutan seed 5 0.985 -> 1.000). (c) **3 clusters good under local mode stop being good** (bar: none): ALLELE-LIKE -> DIVERGED at R 0.996 -> 0.9896, 0.9958 -> 0.9901, 0.9962 -> 0.9893, i.e. across the 0.99042 line; the third has 4 reads, a 4,764 bp consensus trimmed to 934 bp because only that stretch is covered twice; 1 cluster becomes good. Median change of identity x coverage on R +0.00000 (807 clusters), mean +0.00019. The testis E1 cluster is unchanged (30 reads, 5,520 bp, ELSEWHERE, Tm 0.998).
- **T4 not run** (the decision rule is already failed; no seed-9 controls were made).
**Decision: trimming is not added to `discover.py`** (`CONS_MODE = "l"`, no trim). `endtrim.py` stays as an opt-in variant (`lt`). **Reading.** A rule that needs two reads is exact on synthetic tails but, in 3 to 5 read clusters, the window covered twice can be a small part of a real transcript (a long isoform seen once), and it moves clusters across the allele line; the one remaining true tail (ptr_control_s7) is covered twice. The local-mode regression stays at 3 of 745 good clusters (0.4%). Not tried: trimming only the unaligned flank of the consensus (rejected on principle: it would remove exactly the reference-absent exons the discovery path looks for), or reporting an end-support flag per cluster without changing the consensus.

### Amendment 39 (2026-10-08, before any flag is computed): an end-support flag, no change to the consensus

**Why.** Amendment 38 showed that trimming the consensus to >= 2-read support is exact on synthetic tails but moves real clusters across the allele line. A flag leaves the consensus alone and tells the reader which end of a candidate rests on one read.
**Definition (no tuning constant).** Align a cluster's reads (up to 100 longest) to its consensus, primary alignments, `splice:hq -uf`; let `S` be the sorted alignment starts and `E` the sorted ends on the consensus, `n` the number of aligned reads. The 5' end is **single-read** when the extreme read extends past the next one by more than the spread of all the others: `S[1] - S[0] > S[n-1] - S[1]`; the 3' end likewise: `E[n-1] - E[n-2] > E[n-2] - E[0]`. `n >= 3` is required (fewer: flag undefined). Reported per cluster: both flags, the gaps `S[1]-S[0]` and `E[n-1]-E[n-2]` (the single-read flank length), `n`.
**Bars (set now).** (F1, synthetic, sensitivity) the 200 tail seeds of Amendment 38 (301 to 500; 300 bp tail on the longest read, on a random side): **the tailed end is flagged in >= 90%.** (F2, synthetic, specificity) the 200 tail-free seeds (101 to 300): **each end flagged in <= 10%.** (F3, real) of the two class-changing local-mode regressions of Amendment 37 (ptr_control_s7 5-read cluster, ppy_control 5-read cluster), **both are flagged at the tail end** (the tail of the first is covered twice, so this can fail). Reported, no bar: the flagged fraction per class (PRESENT, ALLELE-LIKE, DIVERGED, NOVEL, ELSEWHERE) in the 16 control runs and in the 4 real runs; flank-length distribution; whether the flagged control clusters are enriched among the clusters whose class changes between global and local consensus.
**Decision.** The flag is added to the per-cluster output of `discover.py` iff F1 and F2 hold; F3 is reported either way.

#### Amendment 39 result (2026-10-08; `endtrim.end_flags`, `end_flags_report.py`; `docs/unmapped_rescue_results/end_flags_report.log`, `end_flags.json`): the flag is NOT added to `discover.py`

- **F1 NOT met (sensitivity).** The tailed end is flagged in **29.5%** of the 200 tail seeds (59 of 200; bar 90%); the untailed end in 1 of 200. The 300 bp tail is not an outlier against the natural 5' dispersion of the generator (read starts spread over 500 bp): a read that sticks out by 300 bp where the others spread over ~400 is indistinguishable, from the reads alone, from a longer isoform seen once. The rule works as defined; reference-free, the information is not there at that tail length.
- **F2 met (specificity).** Tail-free seeds: 5' flagged 0 of 200, 3' flagged 2 of 200 (1.0%).
- **F3 NOT met (1 of 2).** ppy_control tail cluster (2,170 bp, 5 reads): flagged at the 3' end, single-read flank 931 bp. ptr_control_s7 cluster (2,650 bp, 5 reads): not flagged (flank 76 bp: its tail is covered twice).
- **Real data, reported.** The flag fires on **21% of PRESENT and 32% of ALLELE-LIKE control clusters** (817 clusters; 25% of all): real read ends have heavy tails, so "the extreme read sticks out" is common and says little about the cluster. Control clusters whose good / not-good status differs between the global and local consensus (69): 46% flagged versus 25% base rate. Flag undefined (< 3 aligned reads): 0.
- **Useful fact for the discovery candidates.** All six real-run clusters with >= 20 reads (human 146, 37, 35, 24 reads; gorilla testis 30 reads, the STON1-GTF2A1L sequence; chimpanzee 24) have a single-read flank of **0 bp at both ends**: the consensus begins and ends where at least two reads begin and end. Of the 106 real-run clusters, 20% of ELSEWHERE and 17% of NOVEL clusters are flagged at either end.
**Decision: not added** (F1 fails). `end_flags` / `end_support` stay in `endtrim.py`; `end_flags_report.py` reproduces the table. Together with Amendment 38 this closes the question for now: the data cannot separate an unsupported tail from a genuine long isoform at the lengths seen, either by trimming or by a flag; what remains available is the descriptive single-read flank length per cluster.

### Amendment 40 (2026-10-08, before the re-run): Amendment 34 on the adopted pipeline

Amendment 34 used the global-mode consensus and the Amendment 33 ELSEWHERE rule. The adopted pipeline is the local-mode consensus (Amendment 37) with rule 35c. **The test, predictions (P1: >= 80% of the gorilla-testis ELSEWHERE reads in same-species clusters; P2: >= 50% of the human ELSEWHERE + NOVEL reads in clusters best on a non-human assembly at >= 0.99) and decision are exactly those of Amendment 34**; only the inputs change: the consensus of every ELSEWHERE, NOVEL and UNSUPPORTED cluster of the four real runs is taken from `cons.l.fa`, classes from `classes.l.a35c.json`. No new bar. The Amendment 34 numbers are superseded where they differ; both are reported.

#### Amendment 40 result (2026-10-08; `SPECIES_VARIANT=l species_check.py`; `docs/unmapped_rescue_results/species_check_local.{log,json}`)

105 consensus sequences (local mode). **P1 met, stronger: 51 of 51 gorilla-testis ELSEWHERE reads (100%) are in same-species clusters** (was 47 of 51); the testis NOVEL clusters (3, 10 reads) still match no gorilla assembly (269 bp: chimp 0.967, orangutan 0.970). **P2 NOT met (0.184 vs bar 0.50; Amendment 34: 0.155).** Human ELSEWHERE 70 clusters / 703 reads (was 56 / 584): CHM13 scores 0.0 on all the top clusters; by reads, 564 (80%) have their best species in the 0.9 to 0.99 band, 113 (16%) are chimpanzee >= 0.99, 23 gorilla, 3 under 0.9. Chimp and orangutan: 0 same-species clusters (6 ELSEWHERE / 50 reads and 11 / 49). The human reading stands unchanged: no CHM13 hit, no panel species at the level of a same-species read.
#### Amendment 40b (2026-10-08, before the run): is "no CHM13 hit" real or a spliced-preset artifact?
**Question.** The human ELSEWHERE / NOVEL consensus sequences are scored on CHM13 with `splice` presets. A non-spliced sensitive preset can align what chaining across introns misses. **Test.** The 77 human ELSEWHERE + NOVEL consensus sequences (local mode) aligned to CHM13 (`target.fa`) with `minimap2 -c -x asm20` (primary records, score = matching bases / consensus length, capped at 1). **Prediction:** >= 50% of the clusters (by reads) reach >= 0.9. **Reading rule (set now):** a median score on the aligned ones near 0.999 would mean human sequence the splice presets failed on (the ELSEWHERE class is then an alignment artifact for those clusters); near 0.98 to 0.99 would mean ape-level divergence from human (a non-human sequence); no hit would mean absent from CHM13.

#### Amendment 40b result (2026-10-08)

**Prediction NOT met: 0 of 77 human ELSEWHERE + NOVEL consensus sequences align to CHM13 with `asm20`** (0 PAF records; the positive control, 6 human control consensus sequences, aligns: a 322 bp exon of a 3,189 bp spliced consensus at 100%, so the run works and asm20 finds exon blocks). Under the reading rule set before the run, this is "absent from CHM13", not a spliced-preset artifact and not human sequence at 0.999 or ape-level 0.98 divergence from it. Together with Amendment 40: sequences that sit at 0.94 to 0.99 of the gorilla, chimpanzee and orangutan assemblies and have no seed in the human assembly, 77 clusters, 739 reads in the human library.
**Source check (descriptive).** All 739 reads come from one SMRT cell, `m64404e_240606_134033`, the same cell as the 905 mapped control reads of that library; chimp and orangutan unmapped reads come from the same cells as their mapped reads; **no cell is shared between libraries**, so this is not read bleed from the other animals' runs. Whether the cell carries more than one sample (multiplexed or mixed-species) cannot be told from the read names. **Unresolved:** the identity of these sequences. A BLAST against nt (external service; not run, the sequences are lab data and sending them needs the owner's agreement) or the sample sheet of the cell would settle it.

### Amendment 41 (2026-10-08, before any truth or pipeline number): recall and false flags on real reads with assembly truth

**Question (user goal).** Take every locus of the maternal assembly that lacks a primary counterpart and has RNA reads; report how many the pipeline flags and how many flags are wrong.
**Setting.** KB3781 = mGorGor1 ("Jim"); RNA = its fibroblast Iso-Seq FLNC aligned to the primary (`fibroblasts/GCA_029281585.2_flnc_mm.bam`, `splice:hq -uf -N 50`). The primary is a per-chromosome mosaic, each chromosome identical to one haplotype (`rna_allele/chrmap.tsv`): primary = paternal on chr 1-4, 7-9, 12-16, 18, 19, 21 (15 chromosomes), maternal on 5, 6, 10, 11, 17, 20, 22, 23, X. Maternal-specific sequence can therefore only lie on the 15 maternal chromosomes paired with a paternal primary chromosome, or on the 201 maternal unplaced / unlocalised contigs (47.2 Mb).
**Truth (assembly-based, read-anchored).** A read is *mat-specific* when (i) its primary alignment on the maternal assembly (`splice:hq -uf`) is clean (de <= 0.00958) and (ii) its *assembly version* (the maternal genomic bases under its aligned blocks, concatenated across introns) scores < 1 - 0.00958 against the whole primary by the pipeline's own score (identity x coverage with terminal-exon rescue, gate off). The read fixes WHICH maternal sequence is expressed; the assemblies decide whether the primary holds it. A *truth locus* = mat-specific reads grouped by overlapping maternal spans (same contig, gap <= 1 kb). Reported with >= 1 and >= 3 mat-specific reads (3 = the pipeline's minimum cluster).
**Candidate reads for the truth (prefilter, loose).** Mask M = maternal bases on the 15 chromosomes with no primary counterpart within delta/2 in the existing primary-vs-partner alignments (`rna_allele/out/chr*.paf`, asm5, cs): bases not covered by any record, plus 500 bp windows whose divergence (mismatches + indel bases) exceeds delta/2; plus all unplaced maternal contigs. P = primary regions hit by M's sequences (minimap2 asm20, -N 50, all records). Candidates = reads unmapped on the primary or whose primary alignment overlaps P. Only candidates are mapped to the maternal assembly and tested. A mat-specific read outside the candidates is missed by the truth (caveat, not measurable here).
**Pipeline (blind to the haplotype assemblies).** Net = every read whose primary alignment is unmapped or has de > 0.00958, genome-wide. Adopted pipeline: seed-round clustering with the proper-overlap edge, star step, abPOA local-mode consensus, support test, identity x coverage on the primary with terminal-exon rescue and the library gate. **Flag** = a supported cluster with R < 1 - 0.00958 (not PRESENT / ALLELE-LIKE). The maternal and paternal assemblies are not used by the pipeline.
**Scoring.** Every read of a flagged cluster is mapped to the maternal and the paternal assembly and given a status: mat-specific, pat-specific (same definition on the paternal assembly), shared (clean on a haplotype and its assembly version has a primary counterpart within delta), or unplaced (clean on neither). A flag is **TRUE** when most of its placed reads are mat- or pat-specific, **WRONG** otherwise. A truth locus is **recovered** when a TRUE flagged cluster has most of its mat-specific reads inside that locus. For every missed truth locus the first failing step is reported: reads not in the net (clean on the primary: absorbed by a paralog or allele within delta), in the net but not clustered, clustered but not flagged, unsupported.
**Positive control (prediction).** The STON1-GTF2A1L readthrough (maternal chr12, about 95.96 to 96.03 Mb; chr12 primary = paternal) is a truth locus with >= 100 mat-specific reads and is recovered.
**Decision (a convention, set now).** "O3 detect + flag works on the real reads of KB3781" is stated iff recall >= 0.5 over truth loci with >= 3 mat-specific reads **and** WRONG <= 10% of flagged clusters. Otherwise the numbers are reported without that claim.

#### Amendment 41b (2026-10-08, after building the mask and the net, before any truth or pipeline number): the truth's candidate reads

**Why.** The registered loose mask (delta/2 per 500 bp window) covers **713.6 Mb** of the 15 maternal chromosomes (456,384 windows of 500 bp with >= 3 divergent bases, mostly heterozygous SNP pairs, plus uncovered segments) and 47.2 Mb of unplaced contigs. Its primary counterpart would contain most of the expressed genome, so "candidates" would be most of the 13.5 M primary reads, which cannot be mapped to the maternal assembly here (time; disk 18 GB).
**Replacement.** Library: 13,516,443 mapped primary reads + 959 unmapped. Net (pipeline input) = 525,813 reads with de > 0.00958 + 959 unmapped = **526,772 reads** (3.9%). Truth candidates = (a) **all 526,772 net reads**, mapped to the maternal assembly and tested exactly as registered; (b) a **1% random sample of the other 12,990,630 reads** (samtools -s 41.01, de <= 0.00958 on the primary), mapped and tested the same way, to **estimate** how many mat-specific reads lie outside the net (reads the net cannot see by construction: clean on the primary, absorbed by a paralog or an allele within delta in gap-compressed divergence). Reads on maternal chromosomes whose primary chromosome IS the maternal one (identical sequence) are mat-specific never; they are not tested.
**Consequences, stated now.** Truth loci are built from the net's mat-specific reads; recall is reported over them ("net-reachable truth"). The outside-net share is reported as an estimate with its count, and the recall over all truth is bounded with it (lower bound: as if every outside-net truth locus were a separate missed locus). The decision convention of Amendment 41 is applied to the net-reachable recall and states the bound beside it. Bars, flag, scoring, positive control: unchanged.

#### Amendment 41 result (2026-10-09; `mat_truth.py`, `run_net.py`, `test_mat_truth.py`; work `/mnt/linuxdisk/tmp/o3_rescue/mattruth/`): the convention is NOT met

**Truth.** Library 13,516,443 mapped primary reads + 959 unmapped; net 526,772 reads. Of the net, 525,921 map to the maternal assembly, **41,121 cleanly (7.8%)**, 41,098 on maternal-specific-capable contigs; **5,331 mat-specific reads** (assembly version < 0.99042 on the primary; 3,825 of them in [0.98, 0.99042), 561 under 0.9, 134 with no primary record). **Truth loci: 940 with >= 1 mat-specific read, 229 with >= 3 (4,511 reads).** Positive control: STON1-GTF2A1L = CM054594.2:95,957,393-96,028,773, **130 mat-specific reads** (prediction met). Outside the net (1% sample, 129,357 clean reads): 436 mat-specific reads (~43,000 genome-wide by extrapolation, ~8 times the net's 5,331); 342 of them fall inside a net truth locus, 94 in **82 further loci** the net cannot see.
**Pipeline (blind, registered settings).** Seed rounds (6 x 2,000 seeds) put 87,848 of the 526,772 reads (16.7%) in 670 components of >= 3; **one component holds 55,967 reads, and the 54 components over 60 reads (kept whole by the star step) hold 80,939**. After the star step: 655 clusters, 86,736 reads. Classes: PRESENT 0, ALLELE-LIKE 69, DIVERGED 518, NOVEL 66, UNSUPPORTED 2. **Flagged: 584 clusters (80,737 reads).**
**Scoring (registered).** **TRUE 102, WRONG 482 (82.5%).** By class: DIVERGED 76 TRUE / 442 WRONG; NOVEL 26 / 40. Of the WRONG flags, 448 have most reads *unplaced* (clean on neither haplotype assembly of the same animal) and 34 most reads *shared*; across all WRONG flags the read mix is 89.6% unplaced, 8.5% shared, 1.9% mat. Even the TRUE flags are 78.8% unplaced reads by count (dominated by the 55,967-read component, TRUE on its 56 placed reads only: 53 mat; its consensus scores 0.619 on the primary). The next two largest flags are real: 6,467 reads (pat 6,451, R 0.973) and 5,249 reads (pat 5,239, R 0.930), paternal-specific sequence on maternal-primary chromosomes.
**Recall.** **35 of 229 truth loci with >= 3 mat-specific reads (15.3%)**; 35 of 940 with >= 1 (3.7%). First failing step of the misses (>= 3): not clustered 188, clustered but held by the primary 4, flagged in a WRONG / other-locus cluster 2. Bound with the outside-net loci of the sample: 35 / (229 + 82) = 11.3% (each sampled outside read stands for ~100, so those loci are taken as >= 3). **STON1-GTF2A1L positive control: recovered.**
**Decision.** Recall 0.153 < 0.5 and WRONG 82.5% > 10%: **"O3 detect + flag works on the real reads of KB3781" is NOT stated.**
**Reading.** Two failures, both at genome scale and both outside what the earlier per-family runs exercised: (1) the seed-round clustering covers 16.7% of the net and chains loci into giant components; 97% of the missed truth loci are never clustered; (2) 92% of the net is reads that are clean on neither haplotype assembly of the animal itself; their clusters have a consensus supported by the reads (UNSUPPORTED 2 of 655) yet about 2% from the primary, and they make most of the flags.
#### Amendment 41c (2026-10-09, before the run): are the WRONG flags noise, or sequence the assemblies lack?
Descriptive, no bar. The 584 flagged consensus sequences aligned to the maternal and paternal assemblies (`splice:hq -uf`, best record, identity x coverage with rescue and the library gate). Reported for TRUE and WRONG flags: the share whose consensus reaches >= 0.99042 on either haplotype assembly (a real haplotype sequence whose reads are noisy) versus not (the consensus itself is absent from both assemblies of the animal: systematic read features, assembly gaps, or another source).

#### Amendment 41c result (2026-10-09; `diag41c.py`; `docs/unmapped_rescue_results/mattruth_{eval,loci}.json`)

| flag verdict | consensus >= 0.99042 on a haplotype assembly | 0.9 to 0.99042 | < 0.9 |
|---|---|---|---|
| TRUE (102) | **97** | 2 | 3 |
| WRONG (482) | 27 | **412** | 43 |

**The TRUE flags are real haplotype sequence** (97 of 102 consensus sequences on the maternal or paternal assembly). **The WRONG flags are not sequence the assemblies lack**: 455 of 482 consensus sequences are on neither haplotype assembly of the animal, and they score the same on the haplotypes as on the primary (median 0.9785 both). Their reads agree with each other (the support test passes) and disagree with every assembly by about 2% per base: a systematic read-level deviation shared by the reads of a cluster (median 6 reads), not a missing copy. 27 WRONG flags have a haplotype-level consensus but noisy reads (973 of the 9,746 WRONG-flag reads).
**What this means for O3.** On real reads the pipeline cannot tell, without the assemblies, a cluster of systematically deviating reads (about 2% from every assembly) from a divergent reference-absent copy (also about 2% from the primary). In KB3781 the first outnumber the second about 4.5 to 1 among the flags, and the seed-round clustering reaches 15% of the expressed maternal-specific loci. The detection step works where it reaches (35 recovered, STON1-GTF2A1L among them) and the classification has no specificity on real reads at genome scale.

### Amendment 42 (2026-10-09, before any code runs): locus-binned clustering at genome scale, same truth, same scoring

**Why.** Amendment 41: the seed rounds clustered 16.7% of the 526,772-read net and chained loci into a 55,967-read component; 188 of 194 missed truth loci were never clustered. A whole-net all-vs-all is 277 billion pairs; the net's reads sit in ~14,700 loci, so clustering inside loci is ~66 million pairs with a 500-read cap.
**Clustering (replaces the seed rounds only).** Each net read is binned by its OWN primary alignment on the primary assembly (contig, strand, overlapping reference spans); the 959 unmapped reads form one more bin. Nothing from the haplotype assemblies is used. Inside a bin: the adopted star step (Amendment 24: proper-overlap edge with every gap counted, delta = 0.00958, longest-first star, clusters of >= 3) on the all-vs-all of the bin's 500 longest unassigned reads; when a bin has more, every other unassigned read is mapped to the star centres just made and joins the first (longest) centre it passes the same edge with; the leftovers repeat the procedure, at most 20 passes. Downstream unchanged: abPOA local-mode consensus, identity x coverage on the primary with rescue and the library gate, support test, flag = supported and R < 0.99042.
**Scoring.** Exactly Amendment 41: the same 940 / 229 truth loci, the same read statuses (maternal and paternal placement of every flagged read, assembly versions on the primary), TRUE / WRONG verdicts, recall, first failing step, STON1-GTF2A1L.
**Predictions.** (P1) recall over truth loci with >= 3 mat-specific reads >= 0.5 (it was 0.153). (P2) the WRONG count rises above 482 and the WRONG share stays above 10%: binning does not touch the systematic-deviation clusters of Amendment 41c. (P3) STON1-GTF2A1L recovered.
**Decision.** The binned clustering replaces the seed rounds for genome-scale runs iff P1 holds and every truth locus recovered by the seed rounds (35) is also recovered or its miss explained. The O3 convention of Amendment 41 is applied unchanged (P2 predicts it fails).

#### Amendment 42b (2026-10-09, mid-run, before any flag or recall number of the binned run): a smaller all-vs-all cap for the 132 largest bins, for compute only

**Seen so far (disclosed).** Bins: 10,331 with >= 3 reads (516,475 reads; largest 61,120; the 959 unmapped reads are one bin). 10,199 bins (299,224 reads) finished with the registered cap of 500: **5.3% of their reads clustered** (the Amendment 24 edge counts every gap, and most net reads are too noisy to pass it with each other), but **45.9% of the truth reads in them clustered**, and 87 of the 183 truth loci (>= 3 reads) lying wholly in finished bins have >= 3 reads in clusters. No consensus, flag, verdict or recall has been computed.
**Problem.** The all-vs-all with base-level alignment costs ~20 s per 400 reads; the 132 remaining bins (217,251 reads; 14 truth loci with >= 3 reads touch them, including STON1-GTF2A1L in the unmapped bin) would take ~16 h at cap 500 with up to 20 passes.
**Change.** For the 132 remaining bins only: cap 200 (all-vs-all of the 200 longest unassigned reads per pass), everything else as registered (centre assignment of the other reads with the same edge, up to 20 passes). The finished bins are kept as they are. A bin that has more than 200 reads now goes through centre assignment where cap 500 would have run the all-vs-all; this can only lose pairs between two non-centre reads. Predictions, scoring and decision of Amendment 42 unchanged.

#### Amendment 42 result (2026-10-09; `binclust.py`, `run_binned.py`, `mat_truth.py` with `MATTRUTH_TAG=_binned`; `docs/unmapped_rescue_results/mattruth_eval_binned.json`)

**Clustering.** 10,331 bins with >= 3 reads (516,475 reads; largest 61,120; the 959 unmapped reads one bin); 10,199 bins at cap 500, the 132 largest at cap 200 (Amendment 42b). **2,894 clusters, 79,897 reads (15.2% of the net)**; no cross-locus chaining (the seed rounds made one 55,967-read component). The Amendment 24 edge counts every gap, so most net reads, which are noisy, never pass it with another read; the truth reads (clean on the maternal assembly) do far more often.
**Classes.** PRESENT 1, ALLELE-LIKE 220, DIVERGED 2,547, NOVEL 124, UNSUPPORTED 2. **Flagged 2,671 clusters (45,803 reads).**
**Scoring (Amendment 41 exactly).** **TRUE 321, WRONG 2,350 (88.0%).** **Recall over truth loci with >= 3 mat-specific reads: 106 of 229 = 0.463** (seed rounds: 0.153); >= 1 read: 110 of 940 = 0.117. Misses (>= 3): not clustered 113, clustered but held by the primary 5, flagged in a WRONG / other-locus cluster 5. STON1-GTF2A1L recovered. Seed rounds versus binned (>= 3): both 32, seed-only 3, binned-only 74; the 3 seed-only loci (8, 9, 23 reads) are in ALLELE-LIKE clusters or not clustered in the binned run.
**Consensus check (as Amendment 41c, descriptive).** TRUE flags: 307 of 321 consensus sequences on a haplotype assembly at >= 0.99042. WRONG flags: 2,238 of 2,350 on neither haplotype (2,157 at 0.9 to 0.99042): the same systematic read-level deviation as in Amendment 41c.
**Predictions.** P1 NOT met (0.463 < 0.5). P2 met (WRONG 2,350 > 482; share 88.0% > 10%). P3 met. **Decision: by the registered rule the binned clustering does not replace the seed rounds** (P1 fails), although it recovers three times as many truth loci with no loss except 3 explained loci. **The O3 convention is not met** (recall 0.463 < 0.5, WRONG 88% > 10%).
**Reading.** Clustering at genome scale is possible on this machine (about 1.5 h of calls for 526,772 reads) and lifts real-read recall to about half of the expressed maternal-specific loci with >= 3 reads; what stays unsolved is specificity: for every real haplotype-specific cluster the pipeline flags about seven clusters of reads that deviate systematically by ~2% from every assembly of the animal, and without the assemblies the two look the same.

### Amendment 43 (2026-10-09, before any feature is computed): a specificity rule for the ~2% systematic-deviation flags, developed on half the genome, tested on the other half

**Problem (Amendments 41c, 42).** Of the 2,671 flags of the binned run, 2,350 are WRONG; 2,238 of them have a consensus on neither haplotype assembly of the animal, about 2% from every assembly, with reads that agree with each other. The TRUE flags (321) are also about 2% from the primary. A rule must separate them using only the reads, the consensus and the primary assembly.
**Split (fixed now).** A flagged cluster belongs to the chromosome of its reads' bin (the primary contig of the bin; chromosome numbers from `rna_allele/chrmap.tsv`). **Development set: odd chromosomes (1, 3, ..., 23) and unplaced contigs. Held-out set: even chromosomes (2, ..., 22), X, Y, chrM and the unmapped-read bin.** All rule design, feature inspection and threshold choice use the development set only; the held-out set is scored once, with the rule frozen and committed first.
**Candidate signals to examine on development (not exhaustive, written before looking):** what the consensus differences to the primary are (substitutions vs indels; indels inside homopolymer runs, a known Iso-Seq error mode; A>G / T>C substitutions, the RNA-editing signature; position near the read ends; unaligned ends vs interior differences), and read-level properties of the cluster (read divergence among themselves and to the primary).
**Bars on the held-out set (set now).** (S1) the rule removes >= 80% of the WRONG flags; (S2) it keeps >= 90% of the TRUE flags; (S3) recall over held-out truth loci with >= 3 mat-specific reads falls by no more than 10% of its value without the rule. **Decision.** The rule is adopted iff S1 to S3 hold. The O3 convention of Amendment 41 is then recomputed on the held-out set with the rule (recall there is bounded by the binned clustering, about 0.46 overall, so the convention can fail on recall alone; that is reported, not hidden).

#### Amendment 43b (2026-10-09, after development-set inspection, before any DNA number): the labels need the animal's DNA

**Seen on the development set (odd chromosomes; held-out untouched).** WRONG and TRUE flags differ from the primary in the same way: substitutions dominate (15.6 vs 9.7 per kb), transition-rich spectra (no A>G excess: 0.36 vs 0.32 of substitutions), homopolymer share of indels equal (0.33), differences spread along the transcript (6 to 7% of clusters with half their differences in one 200 bp window), reads of a cluster not adjacent in the run. WRONG consensus sequences are about 2% from gorilla, human (0.981), chimpanzee (0.981) and orangutan (0.982) alike; TRUE ones are gorilla-like (human 0.973, orangutan 0.944). At 4 of the 8 largest development WRONG loci **every read of the gene** diverges from the primary (median de 0.028 to 0.060): the animal expresses a version of the gene that none of its three assemblies holds. The assembly-based label cannot tell a copy missing from all assemblies (a correct O3 flag) from an RNA-only artifact.
**DNA labels (set now).** For every flagged consensus (all 2,671; the held-out values stay unread until the rule is frozen): the reference segments = its best alignments on the primary, maternal and paternal assemblies (from Amendment 42's alignments) +- 1 kb. **Variant k-mers** = canonical 21-mers of the consensus absent from those segments (`var_all`); `var_pri` = absent from the primary segment only. Each is counted in KB3781's whole-genome reads (SRR26039725, the Amendment 32 tool). **DNA-SUPPORTED** = >= 10 variant k-mers and >= 50% of them seen >= 2 times; **RNA-ONLY** = >= 10 variant k-mers and <= 10% seen >= 2 times; otherwise UNDECIDED (including < 10 variant k-mers). Calibration (positive control): the TRUE flags' `var_pri` k-mers (haplotype-specific, absent from the primary, present in the animal) must be DNA-SUPPORTED in >= 80% of TRUE flags with >= 10 of them; if not, the DNA labels are not used.
**Labels and bars (replace Amendment 43's).** Real = assembly-TRUE, or WRONG with `var_all` DNA-SUPPORTED (genome-present, assembly-absent). Artifact = `var_all` RNA-ONLY. On the held-out set, with the rule frozen: (S1) removes >= 80% of artifact flags; (S2) keeps >= 90% of real flags; (S3) recall as in Amendment 43. Reported: how many WRONG flags are genome-present (development now, held-out after the freeze), the median WGS count of their variant k-mers (lambda_hap = 4.98).

#### Amendment 43b result on the development set and Amendment 43c (2026-10-09, before any held-out number): the frozen rule

**DNA labels (development).** Calibration met: the TRUE flags' haplotype-specific k-mers (absent from the primary) are DNA-SUPPORTED in **164 of 175 (93.7%)** (bar 80%; median 0.895 of them seen >= 2 times in 277.6 M WGS reads). WRONG flags (1,083): **RNA-ONLY 698 (64%)**, UNDECIDED 339, DNA-SUPPORTED 46 (70% of these with a median k-mer count > 25 = 5x a haplotype: repeat k-mers, so their DNA support is doubtful). The ~2% deviation is mostly not in the animal's genome.
**Where the artifacts come from (development).** Best match over a primate panel (gorilla primary, human CHM13, chimpanzee mPanTro3, orangutan mPonPyg2, siamang mSymSyn1 v2.1, the last downloaded from NCBI, md5 checked): RNA-ONLY WRONG flags best on orangutan 279, siamang 244, chimpanzee 144, gorilla 29, human 2. **Most artifact flags are transcripts of other apes in the gorilla library** (cross-sample contamination of multiplexed ape libraries is the obvious candidate; not established here).
**Rules tried on development (registered labels: real = assembly-TRUE or DNA-SUPPORTED; artifact = RNA-ONLY).** R - orangutan > delta: artifacts removed 82.4%, real kept 72.4%. R - max(orangutan, human, chimpanzee) > delta: 98.6% / 26.7%. R > orangutan: 57.9% / 76.9%. Consensus no closer to orangutan than the primary's own sequence (+ delta): 36.7% / 84.6%. Best panel species = gorilla: 95.8% / 64.3%. **Gorilla within delta of the best panel species: 82.2% / 71.9%.** No rule reaches S2 (>= 90% real kept) on development: 50 of the 175 assembly-TRUE flags (gorilla sequence by construction) match human or chimpanzee better than the gorilla primary, the haplotype differing from the primary at that locus (incomplete lineage sorting makes this expected).
**Frozen rule (committed before the held-out run).** Keep a flag iff its consensus' score on the own primary is within delta (0.00958) of the best score over the panel {own primary, human CHM13, chimpanzee, orangutan, siamang} (score = matching bases / consensus length of the best primary-chain record, `splice:hq -uf`). Equivalently: drop a flag that another ape's assembly explains better than the own reference by more than the allele tolerance. **Predictions (held-out):** S1 met (development 82%), S2 NOT met (development 72%), S3 reported. Bars and decision as in Amendment 43b.

#### Amendment 43c result (2026-10-09, held-out scored once; `spec_heldout.py`; `docs/unmapped_rescue_results/{spec_heldout.log, dna_labels.json}`): rule NOT adopted

Held-out flags 1,413 (artifact 838, real 231, undecided 344). **S1 met: 84.4% of artifacts removed** (707 of 838). **S2 NOT met: 80.5% of real flags kept** (186 of 231; assembly-TRUE 131 of 146 = 89.7%). **S3 met: recall over the 114 held-out truth loci with >= 3 reads 0.456 -> 0.421 (relative fall 7.7%).** STON1-GTF2A1L still recovered. Predictions (S1 met, S2 not) held. **Decision: not adopted** (S2).
**O3 convention on the held-out set with the rule:** recall 0.421 (bar 0.5), WRONG share of the kept flags **245 of 376 = 65.2%** (was 1,267 of 1,413 = 89.7%; bar 10%): not met.
**Reading.** The wrong flags of KB3781's real reads are mostly transcripts of other apes (orangutan, siamang, chimpanzee) whose sequence is not in the animal's DNA (variant k-mers absent from 277.6 M WGS reads); a species screen against a primate panel removes ~84% of them and loses 10 to 20% of real flags, because a gorilla haplotype that differs from the gorilla primary can match human or chimpanzee better (incomplete lineage sorting). Specificity on real reads therefore needs either the individual's DNA (decisive here, unavailable for most animals) or a cleaner library; from RNA and reference assemblies alone the precision stays far from the 10% bar.

### Amendment 44 (2026-10-09, before any number): what IsoCon has that the pipeline lacks, and a reconstruction of the part that can bear on specificity

**Diagnosis (IsoCon, Sahlin et al. 2018).** (1) Nearest-neighbour clustering with the best-supported read as centre: partly present (star step, longest read as centre). (2) Iterative correction of reads against their neighbours before clustering: missing; it raises recall (113 held-out misses are never clustered) and would cluster more contaminant reads, so it does not help specificity. (3) A quality-value-weighted significance test of the variants between two candidates: missing (no QVs; Amendment 9's count test); contaminant transcripts are real, high-quality sequence, so it cannot remove them. (4) **Candidates judged against each other, with the read support of each among all reads of the locus**: missing. The pipeline counts only the net's divergent reads and judges a flag only against the reference. IsoCon has no reference-absence or contamination test; nothing in it decides O3 by itself.
**Reconstruction (4).** For every flagged cluster of the binned run with a primary hit: all primary alignments overlapping the hit (at most 300, every n-th read when more; the 959 unmapped reads have no locus) are aligned to two candidates, the flagged consensus and the primary's own sequence under the consensus alignment (its assembly version), `splice:hq -uf`, best record; each read supports the candidate it aligns to (ties and no alignment: neither). `share` = reads on the consensus / reads on either.
**Hypothesis.** A real haplotype-specific sequence of a diploid animal carries about half the reads of its locus; a transcript leaking from another library (Amendment 43: most artifacts are other apes' transcripts) is a small minority. A reference-absent extra copy can be a minority too (it shares the locus with its paralog), so this can cost real flags.
**Test.** Development set (odd chromosomes, unplaced): `share` of real (assembly-TRUE or DNA-SUPPORTED) and artifact (RNA-ONLY) flags, Amendment 43b labels; the rule is chosen there among (a) a one-sided binomial test of `share` against 1/2 (drop when P(X <= reads on consensus | n, 1/2) < 0.01), (b) the same combined with the Amendment 43c species screen (drop when either drops); held-out scored once after the rule is committed. **Bars (Amendment 43b): S1 >= 80% artifacts removed, S2 >= 90% real kept, S3 recall fall <= 10%.** DNA (the animal's WGS, as in Amendment 43b) remains the arbiter; the rule is what can be done without it.

#### Amendment 44 development result and 44b (2026-10-09, before the threshold curve and before any held-out number)

**Development (`locus_support.py`).** Read share of the flag at its locus: **real flags median 0.490 (quartiles 0.278 / 0.667), artifact flags median 0.003 (0.000 / 0.013)**; ~280 assigned reads per locus; every development flag had a locus. Registered rules: (a) binomial minority test (alpha 0.01): artifacts removed 91.8%, real kept 66.1% (assembly-TRUE 131 of 175); (b) (a) AND the species screen: 97.9% / 50.7%. Species screen alone (43c): 82.2% / 71.9%. Neither registered rule reaches S2 on development: the test against exactly 1/2 rejects real alleles with unequal expression at ~280 reads (any share under ~0.43).
**44b: a share threshold, value fixed by a rule set now.** (c) Drop a flag whose share is below tau. **tau = the largest value in {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25} that keeps >= 95% of the development real flags.** No locus: kept. (d) (c) AND the species screen is reported beside it, not chosen. The held-out set is scored once with (c) after tau is committed; bars S1-S3 as in Amendment 43b.

#### Amendment 44b curve and 44c (2026-10-09, after the development curve, before any held-out number)

**44b curve (development, Amendment 43b labels).** tau 0.01 / 0.02 / 0.05 / 0.10 / 0.15 / 0.20 / 0.25: artifacts removed 66.3 / 80.2 / 86.0 / 87.7 / 89.3 / 90.0 / 91.1%, real kept 90.5 / 86.9 / 83.3 / 80.1 / 79.6 / 77.8 / 76.5%. **No tau keeps >= 95% of real flags: rule (c) is undefined as registered.** The real flags with share < 0.02 are 5 of the 175 assembly-TRUE flags (two of them at an 11,898-read locus, every 40th read sampled) and **24 of the 46 DNA-SUPPORTED WRONG flags, whose seen variant k-mers occur 10 to 103 times in the WGS (median per flag), i.e. repeat copy number, not one haplotype (lambda_hap 4.98)**: a k-mer present in a repeat elsewhere cannot show that this sequence is in the genome.
**44c (post hoc to that observation, disclosed; held-out still unseen).** DNA-SUPPORTED additionally requires the median WGS count of the seen variant k-mers <= 5 x lambda_hap = 25 (haplotype-level copy number); otherwise UNDECIDED. Real = assembly-TRUE or refined DNA-SUPPORTED; artifact unchanged (RNA-ONLY). tau is re-chosen by the 44b criterion (largest tau keeping >= 95% of development real flags under the refined labels). Held-out bars S1-S3 unchanged; S2 is also reported for assembly-TRUE flags alone (labels with no DNA ambiguity) and under the 43b labels.

#### Amendment 44c curve and 44d (2026-10-09, before any held-out number): tau frozen

**44c curve (development, refined labels, 189 real flags).** tau 0.01 / 0.02 / 0.05 / 0.10 / 0.15 / 0.20 / 0.25: artifacts removed 66.3 / 80.2 / 86.0 / 87.7 / 89.3 / 90.0 / 91.1%, real kept 94.2 / 92.1 / 87.8 / 84.7 / 84.7 / 83.6 / 82.0%. The 95% criterion again selects nothing (0.01 keeps 94.2%).
**44d (criterion changed on development only, disclosed).** tau = the largest value meeting the held-out bars on development: **tau = 0.02** (artifacts removed 80.2% >= 80%, real kept 92.1% >= 90%; the margin on S1 is 0.2 points, so S1 on the held-out set is a coin toss). **Frozen rule: drop a flag when fewer than 2% of the reads at its locus align better to its consensus than to the primary's own sequence there.** Held-out scored once: S1, S2 (refined labels; also assembly-TRUE alone and 43b labels), S3, the O3 convention, and (d) the rule AND the species screen for comparison.

#### Amendment 44d result (2026-10-09, held-out scored once; `support_heldout.py`; `docs/unmapped_rescue_results/{support_heldout.log, locus_support.jsonl}`): NOT adopted

**Frozen read-share rule (tau = 0.02), held-out:** **S1 met: 80.5% of artifacts removed** (675 of 838); **S2 NOT met: 85.0% of real flags kept** (153 of 180, refined labels; 84.0% under the 43b labels); **assembly-TRUE flags kept 142 of 146 = 97.3%**; the losses are the DNA-supported WRONG flags (11 of 34 kept), sequences the animal's DNA holds at haplotype copy number but no assembly does, i.e. the reference-absent class O3 is after, which sit as minorities at their paralog's locus. **S3 met: recall 0.456 -> 0.447** (relative fall 1.9%). WRONG among the kept flags 286 of 428 = 66.8% (was 89.7%). STON1-GTF2A1L recovered.
**With the species screen as well (reported, not chosen):** artifacts removed 95.0%, real kept 74.4% (assembly-TRUE 87.0%), recall 0.412 (fall 9.6%), WRONG among kept 110 of 237 = 46.4%.
**Reading.** IsoCon has no reference-absence or contamination test; its one idea that bears on specificity, read support among all candidates of a locus, rebuilt here as the share of the locus' reads that a flag carries, separates allele-like real flags (median share 0.49) from leaked transcripts of other libraries (median 0.003). It keeps almost every assembly-confirmed real flag and removes four in five artifacts, but it also removes most genome-present, assembly-absent flags (low-expressed extra copies look like leaks), and the precision stays far from the 10% bar. Missing still: IsoCon's read correction before clustering (recall), and a quality-value variant test (no effect on real contaminant sequence). Only the individual's DNA separates the remaining cases.

### Amendment 45 (2026-10-09, before any code runs): IsoCon-style read correction before clustering, for recall

**Why.** Held-out recall over truth loci with >= 3 mat-specific reads is 0.456 (Amendment 42); 113 of the misses are loci whose reads never pass the clustering edge (every gap counted, NM / block <= 0.00958) with each other because the reads are noisy. IsoCon corrects every read against its nearest neighbours before clustering.
**Correction (fixed now, nothing tuned).** Same bins as Amendment 42 (each read's own primary alignment; the unmapped reads one bin). Inside a bin, reads in length order are cut into chunks of at most 300; one all-vs-all per chunk (`minimap2 -c --cs=long -x map-hifi -N 50 -p 0.1`). Each read is the target of a pileup of its neighbours' alignments (every record with it as target, other reads only); at a column covered by >= 2 neighbours, the read's base, deletion or insertion is replaced when strictly more than half of the covering neighbours carry the same other allele (IsoCon's majority vote; no significance test: with a handful of neighbours the Amendment 9 test cannot fire). One round. Reads of bins of < 3 reads are not corrected.
**Then.** The Amendment 42 clustering (star step, cap 500, and 200 for bins over 500 reads, as run) on the corrected reads; consensus (local mode) from the corrected reads; the support test with the ORIGINAL reads (a consensus must still be supported by raw reads); classification and the Amendment 41 scoring unchanged, truth loci unchanged.
**Bar (held-out = even chromosomes, X, Y, chrM, unmapped bin; the user's bar).** **R1: recall over the held-out truth loci with >= 3 mat-specific reads >= 0.50** (Amendment 42: 0.456). Reported: development-half recall (nothing is tuned, so it is a second test, not a training number), WRONG count and share (expected unchanged or worse: correction also helps leaked transcripts cluster), STON1-GTF2A1L, the share of net reads clustered, and the first failing step of the misses.
**Decision.** Correction is adopted for genome-scale runs iff R1 holds and recall on the development half does not fall.

#### Amendment 45 result (2026-10-09; `readcorr.py`, `test_readcorr.py`, `run_corr.py`, `run_net.py` parallel stages; `docs/unmapped_rescue_results/mattruth_eval_corr.json`): NOT adopted (R1 missed by 2 loci)

**Correction and clustering.** All 516,475 binned reads corrected (chunks of <= 300, majority of >= 2 neighbours; checkpointed per chunk, the 61,120-read bin included). Clustering the corrected reads: **14,495 clusters, 385,373 reads (73.2% of the net**; 15.2% without correction). The consensus and classification stages ran on 4 forked workers (serial and parallel consensus checked identical).
**Classes.** UNSUPPORTED **12,077** (288,271 reads), ALLELE-LIKE 814, DIVERGED 1,464, NOVEL 133, PRESENT 7. The registered support test uses the ORIGINAL reads, and most corrected-read consensus sequences sit more than 0.00958 from their raw reads: correction makes the reads cluster, and the raw-read support test then rejects most clusters. **Flagged 1,597** (2,671 without correction).
**Scoring (Amendment 41).** **TRUE 355, WRONG 1,242 (77.8%**; 88.0% without correction). **Recall over truth loci with >= 3 mat-specific reads: held-out 55 of 114 = 0.482** (without correction 0.456); **development 72 of 115 = 0.626** (0.470); all 127 of 229 = 0.555 (0.463); >= 1 read 134 of 940. Misses (all, >= 3): not clustered 41 (was 113), consensus unsupported 30, clustered but held by the primary 23, flagged in a WRONG / other-locus cluster 8. STON1-GTF2A1L recovered.
**R1 NOT met (0.482 < 0.50; 2 loci short); development recall does not fall (0.470 -> 0.626). Decision: correction is not adopted.** The new loss is the raw-read support test (30 loci); it was registered on purpose and is not changed after the fact on this held-out set.

### Amendment 46 (2026-10-09, before any paternal-truth number): the support test on corrected reads, tested on the paternal truth

**Change.** Amendment 45 unchanged except the support test: the median divergence of the cluster's CORRECTED reads (the reads the consensus was built from) to the consensus, instead of the original reads. Motivation (seen on the maternal truth, both halves): with raw-read support 12,077 of 14,495 clusters were UNSUPPORTED and 30 maternal truth loci were lost at that step.
**Fresh test set.** The maternal truth has been used for recall on both halves; the held-out half three times. **Paternal truth, never used for recall:** the net mapped to the paternal assembly; reads clean there (de <= 0.00958) on paternal contigs that are not identical to a primary chromosome (the partners of the 9 maternal-primary chromosomes 5, 6, 10, 11, 17, 20, 22, 23, X, and the unplaced paternal contigs); assembly version scored on the primary; pat-specific if < 0.99042; loci grouped as before. Recall over paternal truth loci with >= 3 pat-specific reads. Read statuses for the verdicts now come from the maternal and paternal placement of the whole net (one placement for every read), for both variants alike.
**Bars.** (P1) recall over paternal truth loci (>= 3 reads) with corrected-read support >= 0.50; (P2) it is >= the recall of the raw-read-support run (Amendment 45 as run) on the same loci; (P3) the WRONG share of all flags is <= that of the raw-read-support run under the same statuses. **Decision:** read correction with corrected-read support is adopted for genome-scale runs iff P1-P3 hold. Reported: maternal-truth recall of both variants (not a test), STON1-GTF2A1L.

#### Amendment 46 result (2026-10-09; `run_corr_cs.py`, `mat_truth.py eval2`; `docs/unmapped_rescue_results/{eval2_A45_raw_support,eval2_A46_corrected_support,loci.pat}.json`): NOT adopted (P3)

**Paternal truth (fresh).** Net on the paternal assembly: 520,847 mapped, 36,729 clean, 36,692 on testable contigs, **14,783 pat-specific reads; 437 loci, 99 with >= 3 reads**.
**Classes with corrected-read support.** PRESENT 4,954, ALLELE-LIKE 6,981, DIVERGED 2,357, NOVEL 202, UNSUPPORTED 1: most clusters the raw reads had rejected now match the primary within the allele cutoff. Flagged 2,559.
| | raw-read support (A45) | corrected-read support (A46) |
|---|---|---|
| flags; WRONG (whole-net statuses) | 1,597; 1,242 = 77.8% | 2,559; 2,126 = **83.1%** |
| **paternal truth recall, >= 3 reads (fresh)** | 49 / 99 = 0.495 | **56 / 99 = 0.566** |
| paternal, >= 1 read | 59 / 437 | 72 / 437 |
| maternal recall, >= 3 reads (not a test) | 127 / 229 = 0.555 | 148 / 229 = 0.646 |
| STON1-GTF2A1L | recovered | recovered |
Paternal misses (>= 3, corrected support): not clustered 23, held by the primary 16, flagged WRONG / other locus 4 (no unsupported).
**P1 met (0.566 >= 0.50). P2 met (0.566 >= 0.495). P3 NOT met (WRONG 83.1% > 77.8%). Decision: not adopted.** Correction with corrected-read support reaches the recall bar on a truth set never used before, and pays for it in precision: every relaxation that lets real noisy reads through lets leaked transcripts through too.

### Amendment 47 (2026-10-09, before any IsoCon output): does IsoCon itself report the leaked transcripts?

**Question.** IsoCon reconstructs distinct transcripts from RNA alone and keeps a candidate when reads support it significantly. Amendment 43 showed that most wrong flags of KB3781's real reads are other apes' transcripts absent from the animal's DNA. If IsoCon, given a locus' reads as in a targeted experiment, reports those transcripts as supported candidates, then RNA-only reconstruction (IsoCon's published regime) cannot separate them from real copies either.
**Sample (development half, seed 47).** 40 artifact flags of the binned run (Amendment 42): WRONG, RNA-ONLY (Amendment 43b), best panel species not gorilla (Amendment 43 alignments), one per bin; 20 assembly-TRUE flags (controls), one per bin.
**IsoCon input per locus.** All reads of the flagged cluster plus the other primary reads overlapping the flag's primary hit, every n-th to at most 300 in total (the flag's reads always kept); original reads (FASTA, as in the project's earlier IsoCon runs). `IsoCon pipeline` 0.3.3, defaults (min candidate support 2, p 0.01), `--nr_cores 4`.
**Outcome.** IsoCon's final candidates aligned to the flag's consensus and to the primary's own sequence at the locus (`minimap2 -c -x map-hifi`). A locus **reports the flagged transcript** when some final candidate matches the flag consensus at >= 0.99 identity over >= 50% of the shorter sequence and matches it better than the primary's sequence (identity difference > 0.00958).
**Predictions.** (P1) IsoCon reports the leaked transcript at >= 70% of the 40 artifact loci. (P2, control) it reports the haplotype sequence at >= 70% of the 20 TRUE loci. **Reading registered now:** P1 and P2 together = IsoCon keeps real and leaked transcripts alike, so the separation needs DNA or a clean library; P1 failing = IsoCon's test removes leaked transcripts, and its mechanism is to be identified.

#### Amendment 47 result (2026-10-09; `isocon_leak.py`, `iso_leak.sh`; `docs/unmapped_rescue_results/isocon_leak_{loci,scores}.json`)

IsoCon 0.3.3 ran on all 60 loci (median 300 reads; median 16 final candidates at the leak loci, 7.5 at the TRUE loci). **P1 NOT met: IsoCon reports the leaked other-ape transcript as a final candidate at 21 of 40 artifact loci (52%; bar 70%). P2 met: it reports the gorilla haplotype sequence at 15 of 20 TRUE loci (75%).**
**Why P1 fails (descriptive, after the numbers).** The misses follow read support, not origin: at loci where the flagged cluster has more than 4 reads IsoCon reports the leaked transcript at **9 of 12 (75%), the same rate as the TRUE flags (9 of 12)**; at <= 4 reads, 12 of 28 leaked (43%) vs 6 of 8 TRUE. IsoCon's correction and its minimum support act on abundance: a transcript carried by 2 to 4 reads among ~300 is often absorbed; the same would happen to a real, lowly expressed extra copy. When a leaked transcript has the read support of a real one, IsoCon keeps it as often as a real one.
**Reading (as registered, with the mechanism identified).** IsoCon does remove part of the leaked transcripts, by abundance, like the read-share screen of Amendment 44; it has no test of origin, and at equal support it reports foreign and real transcripts at the same rate. RNA-only reconstruction, IsoCon's published regime included, cannot certify that a transcript belongs to the individual's genome in a library that carries other individuals' transcripts.

### Amendment 48 (2026-10-09, user decision, after the Amendment 41-47 numbers): O3 = detection; placement belongs to DNA

**User's framing.** O3's job is to FIND a copy the reference lacks from RNA; whether a found sequence belongs to the individual's genome (and where) is a separate question, answered by DNA. The Amendment 41 convention (recall >= 0.5 AND WRONG <= 10%) mixed the two: its WRONG count is mostly genome membership (Amendments 43-47: transcripts of other apes, absent from the animal's DNA).
**Claim under this framing (all numbers computed before the reframing, so this is a re-reading, not a test).** (1) Erased copies, real reads (Amendments 1-14, held-out bed H, 53 families): 80.1% of the erased copies' reads rescued, 0.07% wrong joins; the consensus of every family cluster lands on its erased copy (40 of 40), at identity x coverage >= 0.999 for 14 of 40. (2) Reference-absent expressed loci at genome scale with assembly truth (correction + corrected-read support, Amendment 46): **recall 0.566 on the never-used paternal truth (56 of 99 loci with >= 3 reads)**, 0.646 maternal; STON1-GTF2A1L found, DNA-confirmed, and found again in a second gorilla. (3) The candidate list is not genome-certified: in KB3781's library about 4 of 5 flags are other apes' transcripts; the read-share screen (Amendment 44) removes 80% of them and keeps 97% of the assembly-confirmed flags; DNA settles the rest.
**What still does not hold under this framing.** Detection of a removed copy is not exact reconstruction (14 of 40 at >= 0.999). Recall on the maternal held-out half stays under 0.5 with raw-read support (0.482); the 0.566 is the corrected-read-support configuration, whose adoption was refused by Amendment 46 on precision, a criterion this framing sets aside. That configuration is the one the claim rests on, and it was evaluated once on a fresh truth set.

#### Amendment 48 deliverable (2026-10-09; `o3_export.py`, `test_o3_export.py`)

The O3 output under the detection framing: `o3_export.py <run_dir> <out_prefix>` writes the flagged candidates of a run as unplaced one-exon transcripts (GTF on their own sequences + FASTA), a BED at the nearest reference locus, and a TSV; attributes: class, supporting reads, identity x coverage on the reference, median read divergence, nearest reference locus, read share at that locus (Amendment 44), `placement "requires_WGS"`. KB3781, configuration of Amendment 46 (`net_run_corr_cs`): **2,559 candidates** (DIVERGED 2,357, NOVEL 202), 2,551 with a nearest locus, 1,526 with read share >= 0.02; gffread extracts all 2,559 transcripts from the FASTA with the GTF. Files in `/mnt/linuxdisk/tmp/o3_rescue/mattruth/deliverable/` (not in the repo: 9 MB). Advisor summary: Claude Docs "O3 Detection Claim".


---

## UNMAPPED_RESCUE_2026-10-08.md

*(original: `docs/archive/2026-10/merged/UNMAPPED_RESCUE_2026-10-08.md`)*

## Rescuing reads that became unmapped: a pool-first, reference-free net (2026-10-08)

> **Independent review (2026-10-08).** Four fresh reviewers audited Amendments 9 to 24. Corrections, limitations and what changed are in the Errata section at the end of `PREREG_unmapped_rescue_2026-10-08.md`; the sentences below that were wrong have been corrected in place.

Goal (user): design a better net or clustering that saves unmapped reads, from a parental haplotype or after erasing a copy, realigning and rescuing the lost reads.
Pre-registration `docs/PREREG_unmapped_rescue_2026-10-08.md` (body + Amendments 1 to 4); code `bench/unmapped_rescue/` (unit-tested); work dir `/mnt/linuxdisk/tmp/o3_rescue/`.
This write-up was rewritten after an independent review that recomputed every number: the held-out numbers reproduced, the dev-bed numbers changed because a BLAST target cap truncated hits (Amendment 4), and several statements were corrected.

### 1. What this is, and what is O3 and what is not

- **O3 (detect and flag a reference-absent copy)** needs a cluster of unmapped reads and its consensus; it does not need a family. Results 4.1 (clustering) and 4.4 (consensus fidelity, real bed) are the O3 results.
- **Family attribution** (sections 4.2, 4.3) puts the reads of a flagged cluster into a known family's net, so that a family-scoped stage can run. It is an O1-style membership question, not O2 (read-to-copy assignment by PSVs), which is not done here.

### 2. The method

1. **Pool first.** Only the unmapped reads, no reference. Each round maps every unclustered read to 2,000 seed reads (`minimap2 -x map-hifi -c`); an edge needs gap-compressed divergence `de` <= 0.00958 (the registered allele cutoff) and an alignment block >= 50% of the shorter read;
   components of >= 3 reads are clusters; the rest go to a new round with new seeds (at most 6 rounds). This replaces the all-vs-all of the prereg body, which did not finish: one shard of 10,400 reads ran past 10 minutes because the 7,000-read families make it quadratic.
2. **Consensus.** One abPOA consensus per cluster (<= 100 longest reads).
3. **Attribute the cluster (optional for O3).** Discontiguous megablast of the consensus against the surviving copies of the family catalog (erased copies removed; `-max_target_seqs 5000`); a family's score is the consensus bases covered by the HSPs of its best copy;
   the cluster goes to the best family iff its score is >= 1.10 x the runner-up's, else it abstains.
4. **Rescue.** Every read of an attributed cluster joins that family's net.

### 3. Beds

| bed | what | pool |
|---|---|---|
| A (dev) | 2026-08-14 whole-genome excision, 162 two-copy families erased from the full genome and the reads realigned | 59,791 deleted-copy reads left unmapped + 959 background reads (plus 1,719 reads without a label, mostly deleted-copy reads the start-only truth rule misses) |
| H (held out for the frozen rule) | A13's 53 held-out multi-copy families, masked genome | 5,313 deleted-copy reads left unmapped + the same 959; no survivor or unlabelled reads at all |
| M (real) | the 959 reads unmapped on the combined primary assembly of KB3781 against its maternal haplotype (124 or 125 map in one 71-kb stretch of mat chr12) | in A's pool, so M's clustering was dev-visible |

### 4. Results

#### 4.1 Clustering (O3; U1 on bed H)

| bed, delta | deleted-copy reads in a cluster of >= 3 | purity | clusters |
|---|---|---|---|
| A, delta/2 | 53,320 of 59,791 (89.2%) | 1.0000 | 240 |
| A, delta (registered) | 57,664 (96.4%) | 1.0000 | 221 (1 to 7 per family, median 2) |
| A, 2 delta | 59,087 (98.8%) | 1.0000 | 203 |
| **H, delta (registered)** | **5,083 of 5,313 (95.7%)** | **1.0000** | 47 (1 to 4 per family); 0 background reads in a family cluster |

Purity is the share of clustered deleted-copy reads whose cluster's majority family is their own. **U1 passes** (bar: coverage >= 0.90 and purity >= 0.99). A has 1 background read in a family cluster.

#### 4.2 Rescue into a family's net with the frozen rule (U2; frozen before bed H was run)

| bed, rule | rescued-correct (share of unmapped deleted-copy reads) | wrong joins (share of joined) | clusters attributed / abstained | cluster accuracy | copies reached |
|---|---|---|---|---|---|
| **H, frozen: cover score, margin 1.10** | 4,255 (80.1%) | 3 (0.07%) | 27 / 20 | 0.963 | 20 of 32 |
| H, strict (margin 1.00) | 4,255 (80.1%) | 474 (10.02%) | 38 / 9 | 0.684 | 20 of 32 |
| H, margin 1.50 | 4,243 (79.9%) | 3 (0.07%) | 25 / 22 | 0.960 | 19 of 32 |
| H, translated mmseqs, strict | 4,257 (80.1%) | 498 (10.47%) | 38 / 9 | 0.658 | 21 of 32 |
| A (dev), frozen | 55,901 (93.5%) | 53 (0.09%) | 126 / 95 | 0.937 | 76 of 118 |
| A (dev), strict | 55,912 (93.5%) | 626 (1.11%) | 179 / 42 | 0.676 | 78 of 118 |
| A (dev), margin 1.50 | 55,861 (93.4%) | 10 (0.02%) | 113 / 108 | 0.973 | 75 of 118 |
| A (dev), translated mmseqs, strict | 46,545 (77.8%) | 9222 (16.54%) | 175 / 46 | 0.589 | 72 of 118 |

- A13's net on bed H joined 557 unmapped reads (10.5%, 99.6% right; register row 1223). Quote the held-out precision as **1 wrong cluster of 27 (cluster accuracy 0.963)**; "3 wrong joins" is that one small cluster, read-weighted.
- The unthresholded top-1 family of the cover score is the true family for 26 of 38 clusters on H (68%; 4,255 of 4,669 reads) and 121 of 185 on A (65%; 55,912 of 56,521 reads): the margin rule removes the 12 and 64 wrong top-1 clusters, which are small.
- **The margin is fragile at 1.10.** Rescued / wrong joins by margin (A then H): A: 1.0: 55,912 / 626 | 1.05: 55,912 / 237 | 1.09: 55,903 / 219 | 1.1: 55,901 / 53 | 1.2: 55,886 / 46 | 1.5: 55,861 / 10 | 2.0: 55,466 / 7; H: 1.0: 4,255 / 474 | 1.05: 4,255 / 138 | 1.09: 4,255 / 129 | 1.1: 4,255 / 3 | 1.2: 4,252 / 3 | 1.5: 4,243 / 3 | 2.0: 4,052 / 3.
  On A the drop from 1.09 to 1.10 is one 166-read cluster (GWFAM316 attributed to GWFAM317 at ratio 1.092); on H it is the 116-read M pile and a 10-read cluster (ratios 1.097 and 1.095). Above 1.10 the wrong joins stay low, but going from 1.10 to 1.50 on A costs 40 reads and 1 copy and removes 43 wrong joins: 1.10 sits just above a handful of wrong clusters, it is not a plateau of identical answers.
- With the cap fixed, the strict rule is wrong at 1.1% of joined reads on A but still 10.0% on H; Amendment 1's account of why a margin is needed was overstated for A.
- Reads of foreign origin joining an attributed cluster are not judged: A has 1,620 (reads without a truth label, mostly deleted-copy reads mislabelled by the start-only truth rule, so A's rescued count is understated by about that many); H has none, so wrong joins for foreign reads are untested there.

**Why 20 to 25% of H's reads stay unrescued** (abstentions by where the true family ranks, H; reads in brackets): no nucleotide hit 8 clusters (414), true family ranked below 2nd 6 (228), true family 2nd 2 (94), true family has no hit while another does 3 (89), 1 background cluster; plus 230 reads in no cluster.
On A: no hit 33 (1,143), ranked below 2nd 43 (283), true family 2nd 1 (157), no hit for the true family 13 (116), true family 1st but margin small 3 (11). The biggest dev cluster that abstained under the capped output (GWFAM133, 3,408 reads) is attributed correctly once the cap is fixed.
The earlier statement that abstentions are mostly near-ties between two close families is wrong: for 9 of the 12 margin or tie abstentions on H the true family is not in the top two.

#### 4.3 Comparators on the same reads (500 seeded unmapped deleted-copy reads per bed; rescued-correct / wrong / abstained)

| bed | B0: map-ont of single reads to the surviving genomic copies, with A13's coverage and divergence floors | B1: single-read dc-megablast, same cover score and margin | pooled, frozen rule |
|---|---|---|---|
| H | 0 / 0 / 500 | 410 / 4 / 86 | 406 / 0 / 94 |
| A (dev) | 0 / 0 / 500 | 475 / 1 / 24 | 469 / 1 / 30 |

- B0 is NOT A13's rule: A13 also aligned the reads to the family's own net reads and joined 557 reads; B0 aligns single reads to the intron-containing genomic copies and finds almost nothing to align (the usable minimap2 records are 33 and 47 reads of 500), so it measures minimap2's inability to align 15 to 30% divergent spliced reads to genomic copies. B1 > B0 is therefore vacuous.
- **U2's first clause passes on H** (rescued-correct 4,255 >= 2,656, wrong joins 3 = 0.07% <= 5%); **its second clause fails**: the pooled method does not beat single-read dc-megablast (406 against 410 on H, 469 against 475 on A).
  So the rescue comes from the sensitive nucleotide score and the margin rule, not from pooling. Pooling is about six times cheaper in wall-clock here (13 minutes against an estimated 80 for per-read alignment, scaled from the sample).

#### 4.4 Consensus fidelity and the real bed (O3)

- **Fidelity (registered metric: identity x coverage >= 0.999 against the unmasked genome).** Over all clusters with a truth family: 14 of 40 on H and 65 of 212 on A; over attributed clusters: 12 of 27 (H) and 39 of 125 (A).
  By identity alone: 20 of 40 and 95 of 212. That every such cluster's best hit overlaps its erased interval (40 of 40, 212 of 212) follows from how the truth family is defined and is not a result.
- **Bed M.** 115 of the 125 reads that map to mat chr12 95.96-96.03 Mb form one cluster (92%; the other 10 are in no cluster) with no other read in it (the held-out bed's pool gives 116); its consensus is 5,694 bp and aligns to the maternal region at identity 1.0000 (5,557 of 5,557 matches, query coverage 0.976).
  **Attribution of this pile is not stable:** against the full 915-copy catalog its best family is GWFAM267 at 407 covered bases against 371 for GWFAM384 (ratio 1.097, just below the margin, so the frozen rule abstains); with the BLAST cap of 500 the same pile was attributed to GWFAM267 (reviewer's recomputation), and against bed A's reduced target set the call was 371 against 369.
  The best hit comes from a short element that hits 831 of the 915 copies (an independent recomputation puts it near 320 bp, 7% of the consensus). There is no family-specific evidence; treat it as unattributed, which for O3 is a complete flag.
- **Run composition.** The 115 reads come from all four SRA runs (8, 16, 82 and 9): not run-exclusive. The dominant run supplies 71% of the pile, close to its 72% share of the 124 chr12 reads that map, and above its 48% share of the 959 background reads.
- **Cost.** Clustering and consensus of the 62,471-read pool took about 10 minutes and the 221 consensus sequences about 2 minutes of alignment; aligning each read alone took 40 s per 500 reads on H and 69 s on A.

### 5. Verdicts against the registered bars and predictions

| bar | verdict |
|---|---|
| U1 (H): coverage >= 0.90, purity >= 0.99 | **PASS** (0.957, 1.0000) |
| U2 (H): rescued-correct >= 2,656 with wrong joins <= 5%, AND pooled > B1 > B0 | **FAIL as a conjunction**: first clause PASS, pooled does not beat B1 (406 vs 410; 469 vs 475 on A), B0 is not A13's rule |
| U3 (M): >= 90% of the 125 reads in one cluster, consensus identity >= 0.99, no other read in it | **PASS** (92%, 1.0000, none), clustering part dev-visible; the family attribution of the pile is unstable (4.4) |
| prediction: U1 holds | held |
| prediction: pooled > B1 > B0, rescued 20 to 60% of N_D | **failed**: pooled is slightly below B1, rescued 80.1% on H |

### 6. Exploratory arms after the held-out run (Amendments 2 and 3; bed H is no longer held out for them)

- **Family union** (the bases covered by any surviving member, a profile proxy): identical to the frozen rule on H and within 3 reads on A.
- **Oracle homology check** of the 13 abstained H clusters with a truth family against only their own family's survivors (blastn word 7, tblastx, permissive e-values): clear homology for 11 (for example 475 bits, E 5e-135, for the 200-read GWFAM158 cluster) against 30 to 150 bits for random families; none for GWFAM244 and one GWFAM348 cluster. The signal exists and the cover score ranks it badly.
- **Repeat-robust specific cover** (a consensus base covered by n families credits 1/n to each; Amendment 2) and with an evidence floor (Amendment 3, f selected on A as the smallest with wrong <= 1%: f = 0.02): A 55,920 rescued, 270 wrong (0.48%); **H 4,325 rescued, 271 wrong (5.90% > 5%): fails on H.** Without the floor H is 457 wrong (9.6%). A floor of 0.05 (not the selected one) gives H 4,321 / 27 wrong and A 55,852 / 14, reported and not claimed.
- **Per-family HMM profile (Amendment 5; HMMER 3.4 from the `hmmer` mamba env; H only, exploratory).** Universe: the 53 H families; profile = MAFFT alignment of the survivors' transcript consensus (148 survivor consensus sequences), `nhmmer` against the 47 cluster consensus sequences, cover score with the frozen margin. In the SAME universe the dc-megablast cover score gives 4,279 rescued / 59 wrong (1.36%) / 21 of 24 copies;
  **the profile gives 4,522 rescued / 125 wrong (2.69%) / 22 of 24 copies** (bit-score variant 4,479 / 183 wrong, 3.9%). By the amendment's condition (more rescued than the restricted baseline with wrong joins <= 5%) this is a gain of 243 reads. Where it comes from: two clusters of families with many survivors, the 200-read GWFAM158 cluster (9 survivors) and the 70-read GWFAM173 cluster (5 survivors), both abstained under the cover score and attributed correctly by the profile.
  What it costs: the GWFAM163 profile attracts small clusters of other families (five clusters of 56, 30, 10, 5 and 3 reads were attributed to it wrongly), 122 wrong reads in nine changed clusters, and the 24-read GWFAM440 cluster that the cover score had right now abstains. Among the 19 families with >= 3 survivors the profile rescues 496 reads with 14.3% wrong joins (baseline 226 with 20.7%), small numbers.
  So a profile helps where a family has many surviving members and needs a specificity guard against hub profiles; two clusters carry the whole gain. It is hypothesis-generating (H is not held out); a new multi-copy erasure bed would be needed to confirm it.

### 6a. Hybrid net (Amendment 6): pool first, then the leftover reads one by one

On bed H, 2,014 of the 6,272 pool reads (32.1%) are not rescued by an attributed cluster; attributing them alone (same cover score, margin 1.10) gives 4,421 rescued (83.2%) against 4,255 for the pooled rule and 4,357 for the B1 estimate scaled to the bed, with 71 wrong joins (1.58% of joined). The extra work is 32.1% of per-read attribution of the whole pool (506 s of about 1,580 s), above the registered 25%, so the registered rule is **not met** and the pool-first net stays "better than A13 (557) and about 6 times cheaper than per-read attribution", not "rescues more than per-read attribution". Bed A was not run: the rule needs both beds, so the H failure already fixes the verdict.

### 6b. The reconstructed loci as extra reference (Amendment 7, bed H)

Consensus sequences from a seeded half of the unmapped deleted-copy reads (plus the 959 background reads) are added to the reference by aligning reads to them and keeping the higher score (the primary a concatenated reference would give). On the other, never-clustered half: **2,505 of 2,657 held-out unmapped reads (94.3%) land on the consensus of their own family at de 0.0011**, absorbed deleted-copy reads move only when the consensus is closer (2,250 of 11,973; 2,247 to their own family; median de 0.1237 -> 0.0031), and **0 of 41,727 survivor reads move**, 0 of 959 background reads move onto a family cluster. All three bars met. This is the augmented-reference test of the earlier NPIPA2 result (row 1223 neighbourhood), now with consensus sequences that were built without any reference.

### 6c. Worked example: LRPAP1 and its gorilla copies (Amendment 8)

With the mother as reference, 381 of 2,875 LRPAP1 reads are poorly placed (de > 0.00958). They form 4 pure clusters (283, 21, 4, 3 reads). The 283-read consensus (1,492 bp) is the paternal c01 (0.9987, the chr12 22.55 Mb copy the mother's assembly lacks within the cutoff); its nearest maternal node is c00 at 0.9886, the registered lift site c01|c03 only 0.9839. The 4-read cluster (2,680 bp) is c01 at 0.9993 with the mother's best 0.9888 and is the one flagged at the registered 0.999 (O3). The 21-read cluster is a c00 allele (0.5% from the mother's), the 3-read cluster matches neither haplotype (0.97). Realigning all reads to the 4 consensus sequences: 84.8% of the held-out net reads move to the consensus of their own copy, and 2.8% of the reads of copies without a consensus move. Details and caveats in the prereg, Amendment 8 result.

### 6d. An IsoCon-style correction step on the consensus (Amendment 9): refuted

Aligning each cluster's reads to its consensus, testing every non-consensus allele against the cluster's error rate (Bonferroni 0.05) and replacing the consensus where a significant allele has the plurality changes the consensus in 53 of 221 clusters on dev A and 10 of 47 on H. Held-out H: identity x coverage >= 0.999 in **15 of 40 clusters against 14 before** (bar: at least 19), and 3 clusters get worse (bar: at most 2). **The decision rule is not met; the consensus-error explanation of the fidelity shortfall is refuted.** Most corrections are exon-scale deletions that trim a mixed-isoform consensus to its plurality isoform, not base repair. The residual edit columns are not allelic either (1.2% and 1.9% are significant minority-variant columns), and the registered metric turns out to be dominated by unaligned consensus ends: on dev A identity alone reaches 0.999 in 95 of 212 clusters, identity x coverage in 65. In the LRPAP1 example the 283-read c01 consensus has identity 1.0000 to pat c01 (0 edits over 1,490 aligned bases); the 0.9987 is 2 unaligned terminal bases, so the flag it missed was a coverage artifact of the registered metric, and the polish changes nothing there.

### 6e. A supervised synthetic world: erase one copy, recover it (Amendment 10)

30 synthetic families (random genes, 5 to 10 exons), three copies each: A and B present, E erased from the reference (divergence of E from its nearest survivor 0.5%, 1%, 2%, 4% or 8%, six families each), 80 reads per copy, exact truth. Recovery works wherever the erased copy is outside the allele cutoff: every E copy from 1% up has a pure cluster, its consensus is at least 0.999 identical to the true transcript (core identity, end gaps excluded) in all 26 clusters that exist (4 of the 0.5% copies have none), attribution is 100% right with 0 wrong, and the erased copy's reads land on its own consensus in 100% of cases. At 0.5% (inside the cutoff) only 12% of the erased copy's reads are in the net, as designed. **One registered bar fails: reads of the surviving copies move onto the erased copy's consensus in 30 to 50% of cases when it is within 2% of a survivor** (raw alignment scores of a spliced genome alignment and an unspliced transcript alignment are not comparable); a divergence rule (move only if the divergence drops) removes it (0 of 960 in every class) and leaves bed H unchanged (Amendment 11, exploratory). The ends: the simulated consensus is never longer than the true transcript and the aligner clips nothing, so jitter, 5' truncation and a polyA tail (which kills the metric entirely, 0 of 26) do not reproduce the real beds; the real beds' unaligned consensus ends are an untemplated 5' run of G (1 to 2 bp in 100 of 212 dev-A consensus sequences, none matching the genome next to it).

### 6f. Divergence move rule adopted; the untemplated 5' G run (Amendments 12 to 14)

The divergence rule is the default of the augmented-reference step: on a fresh synthetic world (seed 20261009, 32 families, three divergence classes not seen before) survivor reads moved 0 of 640 in every class (score rule: up to 320), erased-copy reads landed on their own consensus 100% wherever a pure cluster existed, and bed H and LRPAP1 are unchanged or better (LRPAP1 held-out net reads moving 172 of 191, other-copy reads 3 of 532). The detection edge in the synthetic world sits at the cutoff: 0% of an erased copy's reads leave it at 0.5%, 55% at 0.75%, 75% at 1%, all from 1.5%.
The 5' G run: simulated from the measured distribution (95% of reads, length 1 to 6, mode 2), it reproduces the unaligned 5' clip of the real consensus sequences (22 of 25 against 61% on dev A), and a trimming rule that measures the prefix against the family's surviving copies (T2) restores the registered metric in 21 of 23 synthetic families with no over-trimming and no false trims. **It is not adopted:** on real families the survivors are too divergent to anchor the 5' end (any alignment for only 44 of 221 consensus sequences), so T2 trims 3 of 212 on dev A and none of the 60 attributed clusters that carry the artifact. The artifact will need a library-level signature (the 5' soft-clipped pure-G run on cleanly aligned reads of surviving copies) rather than a per-family reference.

### 6g. The 5' G run read from the library (Amendment 15): works on synthetic, not adopted on real

73% of the cleanly aligned survivor reads of the real library have a 5' soft clip that is a pure G run of 1 to 3 bases (A, T or C: 80 reads in all): the artifact is a property of the library. A gate that tests this (G over-represented among pure 5' clips, binomial p < 1e-6) is open on the real library and on a synthetic world with the artifact and closed on one without. Trimming the leading G run (1 to 3 bases) of every consensus when the gate is open restores the registered identity x coverage metric on the synthetic world (12 to 22 of 24 erased copies, no harm) and raises it on the real beds (dev A 65 to 82 of 212, H 14 to 18 of 40). **It is not adopted:** it fails two of the three registered real-data bars (recall 92% against 95%; templated sequence lost at most 1 base in 67% against 90%), because real 5' ends are GC-rich: 61 of 197 trimmed consensus sequences had no unaligned G at all (105 had no pure-G clip of 1 to 3; the independent review corrected the accounting: templated bases lost at most 1 in 157 of 197 = 80%, still under the 90% bar), so their leading Gs were aligned to the erased copy and the alignment cannot say which are templated and which are artifact Gs that matched the G-rich flank. Practical consequence: the registered 0.999 flag is penalised by this artifact in about half of the clusters, and either the flag should ignore a leading G run of at most 3 bases when the library gate is open (score side, consensus untouched) or the consensus is trimmed and 1 to 3 templated bases are risked.

### 6h. Per-cluster correction of the 5' G run (Amendment 16, rule T4): best of the options, not adopted

Each read's own leading G run is r = j + a: j the templated G's (shared by the cluster) and a the artifact length (a library distribution, read from the aligned survivors' soft clips). The maximum-likelihood j per cluster, then removing only the consensus' leading run beyond j. On a synthetic world with 0 to 3 templated leading G's: the consensus is within one base of the right 5' end in 24 of 24 clusters (T3, trimming everything: 16 of 24), removes templated bases in 1 of 24 (T3: 7), and restores identity x coverage >= 0.999 in 24 of 24 (T3: 20, untouched: 16); the control without the artifact is untouched. Scored against the MODAL read start it fails the registered exactness bars (58% exact, artifact left in 38%); the independent review showed that convention is wrong when the jitter is one base (half the reads start one base in, so a consensus that begins exactly at the true transcript scores +1; the untouched control scores 21/24 exact). With a start-tolerant measure (the consensus begins at ANY observed read start) T4 is exact in 23 of 24 clusters, removes templated bases in 0 and leaves an artifact G in 1 on this world and on a fresh one (seed 20261301: 23, 0, 1; T3 15, 1, 8), i.e. it meets the synthetic exactness bars. It is still not adopted because the dev-A bar fails (below). On real dev A: 65 to 80 of 212 at 0.999 (bar met), but the registered bar that it removes at least the unaligned clip fails (64 of 100): 23 of the 36 shortfalls are clusters of fewer than 5 reads, where the rule does nothing, and 13 are short by one base.

### 6i. The gate-aware flag metric and the IsoCon-style partition (Amendments 17 to 20)

**Gate-aware identity x coverage (adopted, score side).** When the library gate is open, an unaligned leading G run of at most 3 bases leaves the denominator. Synthetic: 12 to 23 and 16 to 24 of 24 at 0.999, controls identical, 0 false passes; real: dev A 65 to 80, H 14 to 17; LRPAP1: the 283-read c01 consensus is 1.0000 to the paternal c01 and is now flagged (mat best 0.9899), the c00 allele cluster is flagged too (the flag never separated an allele from an absent copy).
**Partition into competing candidates (adopted as an opt-in step).** For a cluster: significant variant columns (the Amendment 9 test), pairs of variants that travel together across reads (hypergeometric, Bonferroni), linked blocks (a lone hotspot column is not one), split on the block that divides the reads most evenly, recurse. On a fresh synthetic world (36 families: control, an erased copy with a second isoform skipping an exon, two sibling erased copies 0.5% apart) it splits 0 of 12 controls, separates all 12 isoform clusters into pure leaves (baseline recovers 12 of 24 transcripts, PART 24 of 24) and has no spurious or redundant candidate; the frozen clustering already separated the sibling pairs in that world, while in two earlier worlds (which used the rejected `map-hifi` and edlib alignments) where it joined 5 and 1 pairs PART recovered 24 of 24 (from 19 and 23); a forced-join test with the adopted alignment by the independent review separated 0.5% sibling pairs in 12 of 12 but only 10 of 12 at 0.2%. Two implementation lessons: `minimap2 map-hifi` soft-clips a long deletion (use `splice:hq -uf` and count N as a deletion), and a scorer that takes the best of two orientations by core identity can score 1.0 on garbage (fixed; the first Amendment 18 to 20 verdicts were partly an artifact of it, see the prereg correction). On the real beds PART splits 72 of 212 (dev A) and 19 of 40 (H) clusters and, taking the best candidate per cluster, raises the gate-aware 0.999 count from 80 to 106 and from 17 to 26; it does not split the 22 clusters with more than 200 edits, so they remain unexplained.

### 6j. The 22 dev-A clusters with more than 200 edits (Amendment 21)

They are 3 to 12 reads each, and the reads are fine one by one (baseline divergence 0.0013, MAPQ 60). The consensus is what is wrong: in 16 of the 22 it does not explain its own reads (median read-to-consensus divergence 0.016 against 0.0015 elsewhere). The reads of such a cluster differ 4 to 10 fold in length and start tens of kb apart in the gene (667, 3,657 and 4,508 bases in one 3-read cluster): the edge rule lets a short read bridge two long reads that do not overlap, the component is a chain of partial or alternative-structure reads, and a global POA consensus of a chain is a blend. A truth-free flag (median read divergence to the consensus above the allele cutoff 0.00958) finds 16 of 22 on dev A and 2 of 3 on H, with 4% and 3% of the other clusters flagged; it marks the consensus unsupported and changes nothing. Six clusters on dev A (one on H) are not explained: their consensus fits their reads but differs from the erased copy by 238 to 614 edits. In reads they are 0.2% of the pool; in clusters 10%, and they cap the cluster-level fidelity count (80 or 106 of 212).

### 6k. Chain-aware clustering (Amendments 22 to 24): adopted, opt-in

The chain was a short read CONTAINED in two long reads that disagree: each containment is a proper overlap, and an alignment of two isoforms that differ by an exon-scale deletion has a tiny gap-compressed divergence, so single linkage joins them. Two changes remove it: the edge must be a proper overlap (a suffix-prefix overlap or a containment, not a shared middle), and, inside each component of at most 60 reads, a longest-first star clustering in which a read joins a cluster only through its representative and the pair divergence counts every gap (NM / block). Synthetic chain world (fresh seed, 12 control and 12 two-isoform families at low depth with fragments): isoform transcripts recovered 24 of 24 against 12 of 24, purity 0.91, control untouched. Real beds (corrected after the independent review found that the first star step used minimap2's default all-vs-all limits, which drop edges in components of 20 or more reads; fixed, tested, and re-run): clusters with more than 200 edits 22 to 13 (dev A) and 3 to 0 (H), unsupported consensus 24 to 14 and 3 to 0, erased-copy clusters at 0.999 38% to 52% (152 of 294) and 43% to 56% (31 of 55), at a cost of 1.4 and 2.7 points of read coverage; 71 of the 126 reads of the 22 dev-A clusters (and 14 of 14 on H) leave the clustering, so part of the improvement is removal. The star step touches only components of at most 60 reads (4 to 8% of the reads); larger clusters, where the mixed-isoform consensus problem of Amendment 9 lives, are untouched. Fresh worlds with the fixed step (depth 16, 22 and 50 reads per copy): isoform transcripts recovered 24 of 24 each, controls untouched, purity 0.89 to 0.93. The proper-overlap edge alone is not enough (Amendment 22 failed its bars) and neither is the star step with the gap-compressed divergence (Amendment 23).

### 6l. Specificity controls (Amendment 25)

**PART guard.** A hotspot (a recurrent sequencing error at fixed columns in a quarter of the reads) made PART split 10 to 12 of 12 single-transcript controls; a block-purity test (reads must carry all or none of a block's variants, and the non-carriers must not carry them above the background error rate) brings that to 1, 0 and 0 of 12 on two fresh seeds, keeps exon-skip and 0.5% sibling splits at 12 of 12, and leaves the fresh partition world unchanged (control 0 of 12, isoforms 24 of 24). On the real beds the many-leaf splits disappear (at most 7 leaves against 14) and the split share falls from 33% to 29% (dev A) and 43% to 36% (H), which is still not evidence of distinct transcripts: alleles are not excluded.
**Allele vs copy in the O3 call.** The registered flag (reconstructed on the truth genome, not identical on the reference) flags an allele of a present copy. Requiring a COPY to be more than the allele cutoff (0.958%) from every reference locus separates them: in a synthetic world 12 NULL families give no cluster, alleles at 0.25% and 0.5% give no COPY call (one ALLELE), erased copies are COPY 12 of 12, and on LRPAP1 the c00 allele is ALLELE while both c01 clusters are COPY (the c01 margin is 0.0005). Alleles 1% or more apart are called COPY (4 of 6 at 1%, 6 of 6 at 2%): RNA alone cannot separate them from paralogs.

### 6m. The three-way O3 call on real haplotypes (Amendments 26 and 27)

With the mother's assembly as reference and the father's as truth, the net of the 34 families and LRPAP1 (2,888 reads) gives 37 clusters. The sequence-only call (COPY = reconstructed on the truth haplotype and more than 0.958% from every reference locus) is COPY for GWFAM175_B0 (269 reads, 0.93 to the mother's nearest locus) and for LRPAP1 c01, as predicted, but also for 5 of 28 clusters the study labels 'shared'; with the father as reference it is COPY for 7 of 24 'shared' clusters. These are diverged orthologs: reference and truth hits sit at the corresponding place of the same chromosome, 1 to 7% apart, because real haplotypes of this animal differ that much at a minority of expressed loci. A one-to-one orthology test (map the reference locus' transcript back to the truth haplotype; if it lands on the consensus' truth locus the two are orthologs) turns 9 of those 12 into ALLELE and keeps GWFAM175_B0 and both c01 clusters as COPY; the 3 that stay COPY include a 480-read cluster that is the mirror image of GWFAM175_B0: the two haplotypes each carry a paralog of that tandem array the other lacks (a copy-number polymorphism the study's synteny truth list does not contain). The registered bars were not met in the father's direction (4 of 7, 3 residual calls), so the rule is not adopted; what the data say is that the O3 call needs orthology, not distance. On the erasure bed H, 31 of 55 erased-copy clusters are COPY against the masked genome (the 31 that reach 0.999), none is ALLELE (by construction), 24 are not reconstructed.

### 6n. Discovery from the transcriptome with no truth built in (Amendments 28 to 31)

The 959 reads that do not align to the combined primary assembly (KB3781, real) cluster into three groups of 3 or more reads (120 reads). The 113-read group reconstructs a 5,694 bp sequence that is absent from the primary and identical (1.0000) to a 65-kb locus on the mother's chr12 (95.96 to 96.03 Mb): 84% of the reads that map on either of the animal's other assemblies end up in it. Two further supported groups (4 and 3 reads; 2,445 and 2,330 bp) have no alignment on the primary, the mother's or the father's assembly: candidates for DNA verification (consensus k-mers in the individual's genomic reads); the 839 other reads stay unclustered. Getting there needed two repairs to the classification that the controls exposed: a spliced aligner leaves small terminal exons unaligned (137 bases, two exons of the same locus), so unaligned consensus ends are re-mapped inside the locus window, and a consensus that its own reads do not support is never called. On reads that are clean on both haplotypes the same pipeline calls nothing absent or novel in one control (46 clusters) and, in a fresh one (54 clusters), one 9-read cluster CONFIRMED: the father's version of a locus where the primary carries the mother's, 1.1% apart. That is the limit already measured: an allele or ortholog just over 1% away cannot be told from a paralog by sequence. Candidate sequences: `docs/unmapped_rescue_results/discovery_candidates.fa`.

### 6o. DNA verification of the discovery candidates (Amendment 32)

A 21-mer lookup of the reconstructed sequences in the animal's own whole-genome reads (277.6 M reads, 41.9 Gb, one streamed pass of 4 minutes) separates the discovery run's three clusters. The 113-read, 5.7 kb sequence is in the DNA: 99.8% of its bases lie in 21-mers seen at least twice, median count 6.0, one haplotype at the measured per-haplotype depth of 4.98 (the primary-assembly control sequence is covered 99.9% at 11, two haplotypes). The two small clusters (4 and 3 reads) are not: 0.3% and 0.2% of their 21-mers are present, the same as shuffled copies (0.3%, 0.1%); they align to no gorilla assembly and not to the human genome, and are not GC-poor or repetitive, so they are RNA-only products, not genes missing from the assemblies. The registered covered-bases bar for the shuffled negatives (<= 2%) was exceeded by one control (4.8%) because a few chance k-mers each cover 21 bases; on the k-mer fraction the separation is 97% against 0.3%. Net result of the discovery path on the 959 unmapped reads: one verified reference-absent sequence, two rejected.

### 7. What did not go as predicted, and caveats

1. Pooling does not rescue more reads than aligning each read with the same sensitive score; it is cheaper and gives the consensus.
2. The registered translated attribution (mmseqs) was wrong at 10 to 17% of joined reads; the nucleotide cover score with a margin replaced it on the dev bed (Amendment 1). The margin 1.10 is a dev-bed choice that, on the corrected data, sits just above a few wrong clusters on both beds.
3. Both erasure beds are synthetic (a copy erased from the same animal's genome) and use catalog families as the target set; the real bed M has no family truth, only the maternal alignment.
4. Bed H was run once with the frozen rule (apart from a reporting crash that was fixed and re-run with identical numbers); the code that implements the rule was committed after that run (Amendment 4).
5. The unmapped pile of a novel gene (bed M) is clustered and reproduced, not named.

### 8. Not covered

Poorly placed reads (absorbed on a paralog rather than unmapped); the full 23 GB library as one pool; IG/TR screening of flagged piles beyond the run-composition check; independent replication; a multi-copy erasure bed to confirm the HMM profile result (only H has surviving-member profiles); any change to the shipped binary.


---

## PHASED_TOY_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PHASED_TOY_2026-10-12.md`)*

## Phased-data toy — results (2026-10-12)

Prereg `docs/PREREG_phased_toy_2026-10-12.md` (c5431f41, committed before the run). Code `bench/phased_toy/phased_toy.py`
(seed 7, 1.1 s, 18 MB), outputs `/mnt/linuxdisk/tmp/phased_toy/` (`results.json`, `S*.calls.tsv`, haplotype FASTAs, reads).

Substrate verified by the script: 30-kb CHM13 chr20 slice, SNRPB unit 10,217 bp, V = 44 SNVs (5 exonic; 0.43%), mRNA 1,008 bp.

| scenario | truth | collapsed ref (hap1 only) | phased (hap1 + hap2, mapped separately) |
|---|---|---|---|
| S0 | homozygous | 0 variable columns | 100 SHARED |
| SA | heterozygous allele | 5 variable columns (= the exonic V), alt fraction 0.50 | 100 ALLELE (loci on both haplotypes, partner found) |
| SB | extra copy on hap2 | the SAME 5 columns, alt fraction 0.53 | 50 SHARED (copy A) + 50 EXTRA_COPY (hap2 locus, no partner on hap1) |

P0 0.00 non-SHARED (bar <= .02) PASS · P1 same columns, |diff| 0.03 (<= .10) PASS · P2 SA 1.00 ALLELE / 0 EXTRA, SB 1.00
EXTRA_COPY / 0 ALLELE PASS · P3 150/150 reads on their true haplotype PASS.

SB's 0.53 (not 0.50): three copy-A reads carry a simulated sequencing error at a V column (replayed: the injected
substitutions land on V positions 341/603/676); sampling noise, not copy-B reads.

**Reading.** The collapsed reference receives the same evidence in SA and SB — it can report "a second sequence at 5
columns" but not whether it is an allele or a copy. With both haplotype assemblies, a truth-free rule (which haplotype
explains the read better; does that locus sit between the same flanks on the other haplotype) separates them
perfectly. The toy is the easy case by design: dispersed duplication with unique flanks, 0.43% divergence, 0.1% read
error, perfect assemblies, equal expression, K = 2. Not tested: tandem duplications, gene conversion, K >= 3,
assembly errors, low divergence (PSVs < reads' error budget).


---

## PREREG_phased_toy_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_phased_toy_2026-10-12.md`)*

## PREREG — phased-data toy (2026-10-12)

**Question (collaborator, via the user):** would phased data help? Toy proof of principle only — one real gene on a
very small real chromosome; no claim about real gorilla data.

**Substrate.** CHM13 chr20:2,481,001-2,511,000 (30 kb, real sequence) containing SNRPB (CAT transcript
CHM13_T0135794, 7 exons, minus strand). Copy unit = transcript span +-500 bp (~10.2 kb). Variant set V = 5 exonic
SNVs + intronic/pad SNVs at rate 0.005 (~0.5% divergence, between gorilla nucleotide diversity and the KB3781
p99 exonic haplotype divergence .00958). Seed 7.

**Scenarios (two haplotypes each):**
- **S0 control** — hap2 = hap1 (homozygous).
- **SA allele** — hap2's only copy carries V (heterozygous SNVs). Reads: 50 hap1, 50 hap2.
- **SB copy** — hap2 carries an EXTRA copy B = unit + V inserted at chr20:2,506,000 (dispersed, unique spacer);
  hap2's own copy = hap1's. Reads: 25 hap1-A, 25 hap2-A, 50 hap2-B.

The same V is used in SA and SB and the per-sequence read counts match (50 reference-sequence, 50 V-sequence), so a
collapsed reference receives identically distributed reads in SA and SB **by construction** (allele fraction is not
used as a copy-number signal: allele-specific and copy-specific expression make it unreliable in RNA).
Reads = spliced mRNA, 5'/3' ends jittered 0-30 bp, substitutions 0.001, indels 0.0003 (famsim's settings).

**Arms.**
- **Collapsed** — map all reads to hap1 only (`minimap2 -ax splice:hq -uf`). Report the variable columns
  (alt >= 20% of >= 10 covering reads) and the alt-read fraction. This arm has no rule that can call allele vs copy.
- **Phased** — map reads to hap1 and hap2 separately (as with mat/pat). Truth-free rule per read:
  NM1 == NM2 -> SHARED; else home = haplotype with lower NM, home locus = merged primary spans there;
  ALLELE if the home locus +-2 kb aligns to the other haplotype in ONE alignment (`minimap2 -cx asm10`) covering the
  query end to end (qstart <= 100, qend >= qlen - 100) with a target span within 5% of the query span (the same locus
  between the same flanks, no large indel), else EXTRA_COPY. The 2-kb flank exceeds the unit's 500-bp pad, so it
  reaches unique sequence; on real data the flank must clear the whole duplicated unit (a toy limit).
- A read is **informative** if its window covers >= 1 exonic V position (truth, from the simulator).

**Predictions and bars (decided before any run):**
- **P0** S0: <= 2% of reads non-SHARED.
- **P1** collapsed SA vs SB: identical variable-column sets, alt fractions within 0.10 (the information boundary,
  `NEGATIVE_RESULTS_REGISTER.md:767`, made concrete — true by design, reported to show it, not as a finding).
- **P2** phased: SA >= 95% of informative reads ALLELE and <= 2% EXTRA_COPY; SB >= 95% of informative copy-B reads
  EXTRA_COPY and <= 5% of all reads ALLELE.
- **P3** non-SHARED reads in SA+SB assigned to their true haplotype >= 98%.

**Not claimed.** Real data has unequal expression, larger and tandem duplications (where "which copy is the extra
one" is undefined), gene conversion, K >= 3 copies (`NEGATIVE_RESULTS_REGISTER.md:561`), and needs a
haplotype-resolved assembly of the same animal. Code: `bench/phased_toy/phased_toy.py`; outputs
`/mnt/linuxdisk/tmp/phased_toy/`.


---

## PREREG_missing_copy_posctl_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_missing_copy_posctl_2026-10-14.md`)*

## PREREG — `missing_copy_flag` positive control, re-measured with stable seeds (2026-10-14)

**Written before any scored run.** One plumbing run (3 genes, d = 0.02, seed 999, below) was made to check the commands
on this machine; seed 999 is excluded from every number here.

### Why

`REPRODUCE.md` (missing-copy section) quotes the positive control as "40/40 at 2% divergence, all confirmed", with
"re-measure before quoting": it was measured with the out-of-repo `/mnt/linuxdisk/tmp/gw22/o3/simB.py`, whose mutations
were seeded by Python's per-process `hash()` (B2, `bench/sim.py` docstring), so no read of it can be regenerated. The
2026-09-23 outcome table (`docs/archive/2026-09/PREREG_o3_rna_only_2026-09-23.md`, OUTCOME, positive control):

| divergence | fired / 40 | `reference_absent_candidate` | confirmed (`--confirm` the mutated copies) |
|---|---|---|---|
| 2% | 40 / 40 | 40 | 40 / 40 |
| 1% | 15 / 40 | 15 | 15 / 15 |

This pre-registration re-measures both rows with `sim.py missing-copy genomic` (stable seeds) and the current
`missing_copy_flag`, over five seeds instead of one, and adds a negative control the original did not have.

### Substrate (every input from NCBI, on the 2026-10-14 laptop)

- Genome: CHM13 v2.0 chr20 = NC_060944.1 (efetch, renamed `chr20`; 66,210,255 bp = the T2T length).
- Annotation: RefSeq GCF_009914755.1 (annotation date 2025-08-01), chr20 rows; `utilities gff-to-gtf ... chr20` gives
  4,574 transcripts = the record's `chr20_ref.gtf` count (`docs/DATA.md`).
- Tools: minimap2 2.30-r1287, samtools 1.22.1 (`pixi.toml`); `missing_copy_flag` built `--release` from the commit
  this file is committed in.

**Deviations from 2026-09-23, fixed now:**
1. The home-search index is chr20 only (`minimap2 -x splice:hq -d chm13_chr20.splice.mmi`), not whole-genome CHM13
   (laptop disk and memory). It can only hide an `unannotated_paralogue` call whose better home is on another
   chromosome. The consensus is a 1-2% mutant of the host locus itself, so a better home elsewhere would need a paralogue
   closer to the mutant than its own source; this is not expected to occur, and every `unannotated_paralogue` call is
   listed.
2. Gene IDs are RefSeq `gene-…` IDs (the GTF's `gene_id`); the 40 genes differ from 2026-09-23's in any case (unknown
   seed).
3. Five seeds per divergence (1, 2, 3, 4, 5) instead of one; pooled n = 200 planted copies per divergence.

### Commands (per divergence d ∈ {0.02, 0.01}, per seed s ∈ {1..5})

```bash
python -I bench/sim.py missing-copy genomic chr20_ref.gtf chm13_chr20.fa chr20 40 $d sim_d${d}_s${s} $s --threads 4
missing_copy_flag --bam sim_d${d}_s${s}.bam --fasta chm13_chr20.fa --loci chr20_ref.gtf \
    --index chm13_chr20.splice.mmi --gff chm13_chr20.RefSeq.gff --confirm copies=sim_d${d}_s${s}.copies.mmi \
    --threads 4 --out run_d${d}_s${s}
python bench/score.py missing-copy-control --sim sim_d${d}_s${s} --run run_d${d}_s${s}
```

**Negative control (N), d = 0.02 seeds only:** the `template|…` reads of `sim_d0.02_s${s}.fq` alone (the extra copy's
reads removed), mapped with the shipped `sim.MM2` command, then the same `missing_copy_flag` call (same `--confirm`
index) and the same scorer. Every flag default is used (`--m-min 0.10 --delta-min 0.01 --min-reads 10 --min-sub 3`).

### Unit and scoring

The unit is a **planted gene** (a line of `<sim>.genes.txt`), scored on the table row whose `locus` equals that gene ID.
Overlapping loci (nested genes, readthroughs) share the reads and fire with it (plumbing: 8 rows fired for 3 genes);
they are counted in a separate column, never in the denominators.
- **fired**: the row's `status` starts with `fired`. A gene with no row (fewer than `--min-reads` reads) is not fired.
- **RAC**: fired and `verdict == reference_absent_candidate`.
- **confirmed**: RAC and the `copies` confirmation column is 1 (the confirmation rule unchanged:
  `conf_identity >= 0.99` and `conf_identity - host_identity >= delta/2`).

### Bars (pooled over the five seeds)

| id | divergence | measure | bar | 2026-09-23 value |
|---|---|---|---|---|
| **P1a** | 2% | fired / 200 | ≥ 190 (0.95) | 40 / 40 |
| **P1b** | 2% | RAC / fired | ≥ 0.95 | 40 / 40 |
| **P1c** | 2% | confirmed / RAC | ≥ 0.95 | 40 / 40 |
| **P2a** | 1% | fired / 200 | in [50, 100] (0.25-0.50) | 15 / 40 = 0.375 |
| **P2b** | 1% | confirmed / RAC | ≥ 0.95 | 15 / 15 |
| **N** | template reads only (2% seeds' genes) | fired / 200 | ≤ 2 (0.01) | not measured |

Per-seed values are reported next to the pooled ones; the bars apply to the pooled values only.

### Reading, fixed now

- **P1a-c and P2a-b and N pass** → the control reproduces: `REPRODUCE.md` replaces "40/40 … re-measure before quoting" by
  the pooled numbers and this file's commands.
- **P1 fails** → the 40/40 does not reproduce on stable seeds; `REPRODUCE.md` and `docs/THESIS_OBJECTIVES.md` quote the
  new numbers, and the failing step (mixture, consistency, verdict or confirmation) is named from the scorer's funnel.
  No threshold of `missing_copy_flag` is changed in response (that would need its own pre-registration).
- **P2a outside its interval** → reported as a change in the 1% sensitivity edge, with the direction; not a failure of
  the detector by itself (the 2026-09-23 edge rested on one seed of 40).
- **N fails** → the detector fires on error-only reads at a simulated locus; reported as a false-positive rate on
  simulated negatives and investigated before any O3 number from `missing_copy_flag` is quoted.

### Plumbing run (seed 999, excluded)

3 genes at d = 0.02 (gene-CST4, gene-STK35, gene-SYS1-DBNDD2): 80 reads, simulation 5 s, flag about 2 s; the three
planted genes' rows and five overlapping rows fired. This checks the commands only and is no part of any bar.

## OUTCOME (2026-10-14) — every bar passes; the control reproduces on stable seeds

`missing_copy_flag` built `--release` at 49e62be3 (binary sha1 159643eb), `bash bench/o3_posctl.sh REF OUT` under
`pixi run`; md5 `chm13_chr20.fa` 1d0b5621077818120e7d40dae01ad4d7, `chr20_ref.gtf` 84145d8099109ddb1b58afc059017bc2.
Run products: `~/rustle_data/o3_posctl/scored/` on the laptop (one `*.per_gene.tsv` per run).

| seed | 2% fired | RAC | confirmed | N fired | 1% fired | RAC | confirmed |
|---|---|---|---|---|---|---|---|
| 1 | 39 | 39 | 39 | 0 | 13 | 13 | 12 |
| 2 | 40 | 40 | 39 | 0 | 21 | 21 | 20 |
| 3 | 40 | 40 | 40 | 0 | 21 | 21 | 21 |
| 4 | 40 | 39 | 39 | 0 | 19 | 18 | 18 |
| 5 | 40 | 40 | 40 | 0 | 17 | 17 | 17 |
| **pooled** | **199 / 200** | **198** | **197** | **0 / 200** | **91 / 200** | **90** | **88** |

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| P1a | 2% fired / 200 | 199 (0.995) | ≥ 190 | ✅ PASS |
| P1b | 2% RAC / fired | 198 / 199 (0.995) | ≥ 0.95 | ✅ PASS |
| P1c | 2% confirmed / RAC | 197 / 198 (0.995) | ≥ 0.95 | ✅ PASS |
| P2a | 1% fired / 200 | 91 (0.455) | 50-100 | ✅ PASS (above the 2026-09-23 0.375, inside the interval) |
| P2b | 1% confirmed / RAC | 88 / 90 (0.978) | ≥ 0.95 | ✅ PASS |
| N | template-only fired / 200 | 0 | ≤ 2 | ✅ PASS (all 200 `no_mixture`) |

**Reading (as fixed):** the control reproduces. The 2026-09-23 "40/40 at 2%" becomes **199/200 fired, 198
`reference_absent_candidate`, 197 confirmed (five seeds)**; the 1% edge is **91/200 fired, 88/90 confirmed**; the detector
fired on none of the 200 planted genes when the extra copy's reads were removed. No `unannotated_paralogue` call occurred,
so deviation 1 (chr20-only home index) hid nothing on this substrate.

**The misses, named (reported, nothing re-tuned):** every one is at a lncRNA / uncharacterized `LOC` gene.
- *Mixture* (1, at 2%): `gene-LOC105372546` (seed 1), `no_mixture`.
- *Consistency* (2): `gene-LOC107985409` (2%, seed 4: n_psv 2, shared 0.26) and `gene-LOC105372659` (1%, seed 4: n_psv 4,
  shared 0.37), `scattered`.
- *Confirmation* (3): `gene-LINC00659` (2%, seed 2: conf 0.9943 vs host 0.9915, margin below delta/2 = 0.0067),
  `gene-LOC105372727` (1%, seed 1: conf 0.9886 < 0.99, host 1.0000) and `gene-LOC107985409` (1%, seed 2: conf 0.9936 vs
  host 0.9893). A host identity of 1.0000 at delta 0.02 says the patched consensus covers almost none of the PSVs (it is
  built over the template read's span only); the confirmation rule is then right to withhold.


---

## PREREG_o3_candidates_foreign_check_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_candidates_foreign_check_2026-10-14.md`)*

## PREREG — are `candidate_copies`' remaining false flags other apes' transcripts? (2026-10-14; to run on the main machine)

**Written before any of the alignments below.** It follows the F < 5 branch of `docs/PREREG_o3_fewread_flags_2026-10-14.md`:
a clean simulation (96 NULL-type units, paralogues down to 1%, present copies down to 4 reads, 5' G read ends) produced no
few-read false flag, so the A14 / H false flags come from something the simulation lacks. This file tests the candidate with
independent evidence on the same library (M3).

### Background, and what is and is not known

- The A14 (no-deletion control) and H (held-out, nothing deleted) runs use reads of the KB3781 fibroblast library. After the
  Amendment 15 fix they still carry false flags: A14 16 (2 b_allele + 14 c_unmatched, 10 families; C1' bar ≤ 8 failed), H 10
  (2 b + 8 c, 6 families; H1 bar ≤ 4 failed) (`docs/O3_CANDIDATES_AMENDMENT15_2026-10-05.md`).
- On the same library, `docs/PREREG_unmapped_rescue_2026-10-08.md` Amendment 43 found that most RNA-only wrong flags of the
  unmapped-rescue chain best match another ape (orangutan 279, siamang 244, chimpanzee 144, gorilla 29, human 2 of 698), with
  their variant k-mers absent from the animal's WGS. How such reads entered the library was NOT established (cross-sample
  contamination of multiplexed ape libraries is the stated candidate; no SMRT cell is shared with the project's own chimp /
  orangutan libraries). Those flags are not the `candidate_copies` flags; whether the A14 / H false flags are of the same kind
  has never been checked.

### Inputs (main machine; nothing here is re-run)

- Flagged unions: `iso.contigs.fa` (+ `iso_names.tsv`) of the corrected runs, `/mnt/linuxdisk/tmp/rna_allele/a14a15/` and
  `/mnt/linuxdisk/tmp/rna_allele/a15h/` (`accept_o3_candidates.sh contigs`), with each union's Amendment 9 class from the same
  runs' classify output (a / b / c / pri). The unit is a flagged union; **false** = class b, c or pri (16 + 10 = 26 expected;
  any other count is reported and used as found).
- Genomes, each aligned as the FASTA or the splice index already on disk:
  gorilla primary `GGO.fasta` (GCF_029281585.2) and KB3781's maternal / paternal assemblies (`mGorGor1.{mat,pat}.splice.mmi`)
  = **own**; human CHM13 v2.0 (`npip_ladder/idx/target.splice.mmi`), chimpanzee mPanTro3 (GCF_028858775.2), orangutan mPonPyg2
  (GCF_028885625.2), siamang mSymSyn1 v2.1 (GCF_028878055.3, the copy Amendment 43 used) = **foreign**.

### Command and score

`minimap2 -c -x splice:hq -uf -N 20 -t 4 GENOME iso.contigs.fa > iso.<genome>.paf` per genome (minimap2 2.30).
Score of a union on a genome = matching bases of its primary records / union length (capped at 1; Amendment 43c's score).
`own` = max over the three gorilla assemblies; `foreign` = max over the four others, with the species that gives it.

### Classes, fixed now (per false union)

- **FOREIGN**: foreign ≥ 0.995 and foreign − own ≥ delta/2 = 0.00479 (Addendum 1's foreign rule, 2026-09-23).
- **FOREIGN-LEANING**: foreign − own ≥ 0.00479 but foreign < 0.995.
- **OWN**: own ≥ foreign − 0.00479 (no other ape explains it better; the next suspects are real-read error structure and
  real paralogue complexity).

**Self-check, required:** the class-a flags (GWFAM175_B0 in A14 and the H class-a union) must come out OWN with own ≥ 0.999.
If either comes out FOREIGN the score is broken and nothing below is read.

### Reading, fixed now

| FOREIGN among the false unions (A14 + H pooled, n ≈ 26) | reading |
|---|---|
| ≥ half | **M3 is the main cause.** The stage's remaining false flags are other apes' transcripts in the library. A species-panel screen becomes the candidate fix, under its own pre-registration that counts its cost in real flags (Amendment 43c's screen kept only 80.5% of real flags); C1' / H1 recomputed without the FOREIGN unions are reported beside, never as a pass. |
| ≤ 5 | **M3 is not the cause.** The false flags are gorilla sequence; each is diagnosed individually (read origins against mat / pat, per-column cover of its consensus, the haplotype it best matches) before any rule is proposed. |
| in between | mixed; both follow-ups, per union. |

Reported beside: every false union with its length, reads, d, own and foreign scores and best species; FOREIGN-LEANING
counts; the same table for the class-a and true flags.

## OUTCOME

(to be appended after the run)


---

## PREREG_o3_fewread_flags_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_fewread_flags_2026-10-14.md`)*

## PREREG — `candidate_copies` few-read false flags: does a clean simulation produce them, and by which mechanism? (2026-10-14)

**Written before any `candidate_copies` run on the worlds below.** One generation-only probe (`synth_world.py build
--spec-world --seed 7`, no alignment, no stage run) was made to read the world's file layout; seed 7 is used nowhere below.

### The open item

Register 1247 / `docs/O3_CANDIDATES_AMENDMENT15_2026-10-05.md`: after the Amendment 15 consensus fix the opt-in stage
still fails both no-deletion controls — A14 C1' 10/53 families with a b/c/pri flag (bar ≤ 8), held-out H1 6/30 (bar ≤ 4).
The remaining false flags are described as "few-read (median 11.5 reads), mismatch-dominated consensus (median whole-length
d 0.036, median id x cov 0.974) — a future amendment's subject". Two mechanisms are suggested by the code and have never
been measured:
- **M1 mixing**: single linkage at de ≤ delta = 0.00958 chains reads of two copies; the column-wise vote gives a mosaic
  more than delta from every locus.
- **M2 consensus error**: reads of ONE present copy, but columns covered by < 3 members are never voted (`VOTE_MIN_COVER`)
  and ends covered by 2 survive the trim, so a 3-12-read cluster keeps its template read's errors and end artifacts.

A third explanation is not a stage defect: **M3 foreign transcripts**. Both controls use the KB3781 fibroblast library, and
the later unmapped-rescue study (`docs/PREREG_unmapped_rescue_2026-10-08.md`, Amendments 41-46) found about 4 in 5 of
that library's reference-absent flags to be other apes' transcripts. A median whole-length d of 0.036 is in the range of a
gorilla-orangutan transcript difference. Such reads ARE absent from the gorilla reference, so no RNA-internal rule can
separate them from a real absent copy at equal support; only a foreign-genome screen can. A simulation contains no M3 by
construction — which is what makes it the right instrument to ask whether M1/M2 exist at all.

### Phase 1 (this file): a diagnostic, not a rule

**Worlds.** `bench/o3_fewread/world.py` builds, per seed, one synthetic world with `bench/famsim` through
`synth_world.build_world` (random 60-kb background per family, synthetic 5-10-exon genes, the shipped HiFi read model
err 0.001 / indel 0.0003), 48 families:

| type | n | copies (in the reference unless marked) | expression (reads per copy) |
|---|---|---|---|
| NULL3 | 8 | A, B (2% SNPs), C (3%) | 80 each |
| NULLclose | 8 | A, B (1%), C (1.5%) — inside single-linkage chaining reach | 80 each |
| NULLlow | 8 | A, B (2%), C (3%) | A, B 80; C one of 4, 6, 8, 12, 16, 24, 32, 48 (one per family) |
| ALLELE@d | 10 | A, B (3%), H = A + d SNPs, **H absent** (d = 0.25, 0.5, 1, 2, 4%; 2 each) | 80 each |
| ERASED | 14 | A (2%), B (3%), **E absent** | A, B 80; E one of 4, 6, 8, 12, 16, 40, 80 (2 each) |

Read variants (both run): **E0** (famsim reads, 5' jitter 30) and **E3** (jitter 3 + the untemplated 5' G runs measured on
the real reads, Amendment 13: `synth_world.add_g5`). Seeds: **1 and 2** for Phase 1. Seeds 3, 4, 5 are reserved, unbuilt,
for a Phase 2 held-out test.

**Pipeline per world × variant.** Reads mapped with the shipped `sim.MM2` command to `genome.ref.fa`; a copies table
(`P.fam.copies.tsv` / `.fa` in `mcl_families`' layout: one family per synthetic family, one copy per IN-REFERENCE copy,
its truth exons, its spliced exon sum) built from the world's truth; `candidate_copies` at its defaults
(`--delta 0.00958 --min-cluster 3 --min-support 6`); minimap2 2.30, binary built `--release` from this file's commit.

**Scoring (`bench/o3_fewread/score.py`), per flagged candidate:** its reads (`reads.tsv`) mapped to their true copy (read
names carry it); `top_frac` = the share of the most common source copy; that copy's kind; the candidate's union sequence
(`contigs.fa`) aligned to every true transcript (`minimap2 -c -x splice:hq -uf`), best hit by identity x coverage. Class:
- **true**: ERASED family, `top_frac` ≥ 0.9 from E (the planted absent copy).
- **allele**: ALLELE family, `top_frac` ≥ 0.9 from H (correctly detected, but class b for C1'/H1).
- **M1 mixed**: `top_frac` < 0.9.
- **M2 pure-present**: `top_frac` ≥ 0.9 from an IN-REFERENCE copy.

A flag is **few-read** iff `n_reads` ≤ 20 (the real false flags' median is 11.5).

**Reported:** per type × variant × seed, families with ≥ 1 flag by class; every M1/M2 flag with n_reads, n_clusters, d,
top_frac, best truth hit; ERASED recall (families whose E is flagged) by E's expression; the n_reads of every true flag (the
cost curve any read floor would impose).

### Branching rule, fixed now

Let **F** = the number of (family, variant, seed) units of the three NULL types (8 × 3 × 2 variants × 2 seeds = 96 units)
carrying ≥ 1 few-read M1 or M2 flag.

- **F ≥ 5 (≥ 5.2%; the real H1 bar is 4/30 = 13%)** → the simulation reproduces the regime. Phase 2 = a NEW amendment to
  this file, written before seeds 3-5 are built, fixing (a) the rule chosen on seeds 1-2 only, (b) its held-out bars on
  seeds 3-5 (NULL-type false-flag units and ERASED recall per expression level), and (c) the real-data acceptance on the
  registered A14 / H substrates with the unchanged C1' / H1 / A13-1 bars, to run on the main machine.
- **F < 5** → the clean simulation does not reproduce the few-read false flags; no rule is developed on simulation (it
  would be fitted to nothing). The next step is the **M3 cause check on the real flags** (main machine, files
  `/mnt/linuxdisk/tmp/rna_allele/a14a15/`, `a15h/`): align the 16 + 10 false-flag unions to CHM13, mPanTro3 and mPonPyg2
  with the Addendum-1 foreign rule (coverage ≥ 0.8, identity ≥ 0.995 and above the gorilla host identity by ≥ d/2), pre-
  registered separately before that alignment is made.
- Whatever F is, any M1 or M2 flag is diagnosed individually (consensus vs template, per-column cover) and reported; a
  defect found that way gets its own fix and regression test, not a threshold.

Nothing in `candidate_copies` is changed before Phase 2 is registered.

## OUTCOME (2026-10-14) — F = 0: the clean simulation does not reproduce the few-read false flags

`candidate_copies` built `--release` at 49e62be3; `bench/o3_fewread/world.py build|run` and `score.py` under `pixi run`
(minimap2 2.30, samtools 1.22.1). Each world × variant: 48 families, about 130 s, peak RSS 0.24-0.34 GB. Run products:
`~/rustle_data/o3_fewread/w_s{1,2}/{E0,E3}/` on the laptop (`flags.tsv` = one row per flagged candidate).

| type (per world) | E0 s1 | E3 s1 | E0 s2 | E3 s2 |
|---|---|---|---|---|
| NULL3 / NULLclose / NULLlow: families with any flag | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| ALLELE (10): allele flags | 4 | 4 (+1 M1) | 6 | 6 |
| ERASED (14): E flagged | 12 | 12 | 12 | 12 |

- **F = 0 of 96 NULL-type units** (no M1 or M2 flag of any read count in a NULL family, including NULLclose's 1-1.5%
  paralogues inside chaining reach and NULLlow's present copies at 4-48 reads). **Branch: F < 5.** No rule is developed
  on simulation; the next step is the M3 cause check on the real flags (`docs/PREREG_o3_candidates_foreign_check_2026-10-14.md`).
- **The one M1 flag** is not few-read: `f29` (ALLELE@0.01, seed 1, E3), one cluster of 131 reads = 80 of the absent allele H
  + 51 of the present copy A (H is 1% from A, just above delta = 0.00958), best hit H at id x cov 0.9979. A chained
  allele+copy cluster, reported, not a few-read flag.
- **ERASED recall by E's reads (pooled over 4 units):** 4 reads 0/8 (below `--min-support 6`, by design), 6-80 reads 48/48.
  All 48 true flags are read-pure (top_frac 1.0), best hit their own E, id x cov ≥ 0.9969 (32 at ≥ 0.999). The 32 true flags
  at ≤ 20 reads have id x cov ≥ 0.9978: **on clean reads a 6-16-read cluster gives an accurate consensus**, so a read floor
  above 6 would remove real low-expression copies here without removing any false flag.
- **Allele flags:** 20, every one at d ≥ 0.01 (above delta), none at d ≤ 0.005 (linked, as designed); best hit H in 20/20.

**Reading.** On famsim reads (uniform HiFi errors, synthetic genes, one or two paralogues), neither mechanism M1 nor M2
produces few-read false flags. The real ones therefore come from something the simulation lacks: foreign or contaminant
transcripts (M3), real-read error structure (homopolymer indels, chimeras, incomplete processing), or real paralogue
complexity (many copies, gene conversion). M3 is the one with independent evidence on the same library, and it is directly
testable on the real flags; it goes first.


---

## PREREG_o3_structure_projection_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_structure_projection_2026-10-14.md`)*

## PREREG — O3 structure projection: held-out simulation test (2026-10-14)

**Written after development on seeds 11, 12, 13 and before any held-out world is built.** Design and the frozen rule:
`docs/superpowers/specs/2026-10-14-o3-structure-projection-design.md` §3 (reference implementation `bench/o3_structure/proto.py`,
md5 511562b804907e1f3cb54b008000dc68; worlds `bench/o3_structure/world.py`, md5 fd2bc142f13a18ee4acfb63f95d05982).

### Claim under test

`candidate_copies`, given the reads of a family with one expressed copy absent from the reference, outputs a candidate that carries
that copy's **exon structure**: every exon boundary of the true copy inside the candidate's span, and no other, within ±5 bp — for
copies that differ from their nearest reference paralog by SNPs (≤ 4%) or by a structural change (exon loss, gain, duplication,
shuffle, 5' truncation, resized intron), on either strand.

### P0 — implementation parity (gate, before the held-out worlds are built)

The Rust implementation in `candidate_copies` (product `<out>.structure.tsv`) must give the **same boundary list as the frozen
prototype for every flagged candidate** of the dev worlds s11, s12, s13 × E0, E2 (148 candidates). A difference is a defect of the
implementation, fixed in the implementation (never in the rule), and P0 re-run. The held-out worlds are built only after P0 passes.
Also required: `cargo test --release` fails only the known tests (the 3 `pipeline_driver_dispatch` tests and, on this laptop, the 10
missing-fixture tests), and the products that existed before (`candidates.tsv`, `clusters.tsv`, `contigs.fa`, `families.tsv`,
`reads.tsv`, `clusters.fa`, `nets.fa`) are byte-identical to the pre-change binary on the dev worlds.

### Held-out worlds

`world.py build OUT --seed S --strands` for **S = 21, 22, 23, 24, 25**; read variants **E0** and **E2**; `world.py run`
(minimap2 2.30, `candidate_copies` at its defaults). 27 families with an absent copy E per world (DIV@0.5/1/2/4/8% × 2, LOSS × 3,
GAIN × 3, DUP × 2, SHUFFLE × 2, TRUNC5 × 2, INTRON × 2, LOWEXPR 8/12/16 reads) and 4 NULL families.

### Unit and scoring (`proto.py WORLD VARIANT --score-rust`, the prototype's scorer applied to the Rust product)

Unit = a flagged candidate of an E family whose union's best true-transcript hit (`minimap2 -c --eqx -x map-hifi`) is that family's E
on `+` ("scored"). True boundaries = E's exon junctions in transcript order inside the union's aligned span (≥ 20 bp from its ends),
mapped to union coordinates through that alignment. A candidate is **structure-correct** iff every true boundary has a projected
boundary within ±5 bp and every projected boundary is within ±5 bp of a true one. Boundary recall / precision pooled over boundaries.

### Bars (pooled over 5 seeds × 2 variants)

| id | measure | bar | dev value |
|---|---|---|---|
| **S1** | structure-correct / scored, all scenarios except DIV@0.08 | ≥ 0.90 | 0.986 (all scenarios incl. 8%) |
| **S2** | structure-correct / scored, per structural scenario (LOSS, GAIN, DUP, SHUFFLE, TRUNC5, INTRON) | ≥ 0.80 each | 1.00 each except SHUFFLE 5/6 |
| **S3** | boundary recall and boundary precision (±5), all scenarios except DIV@0.08 | ≥ 0.97 each | — (≥ 0.99) |

DIV@0.08 is reported, not barred (8% is beyond the family-edge identity the stage's families are built at; dev 6/6).

### Reported beside (no bar)

Exact (±0) structure-correct rate; every scenario and strand separately; candidates not scored and why; detection (E flagged) per
scenario and expression; flags in NULL families; **reconstruction** — identity x coverage of the union to E (≥ 0.999 = exact), and the
same for **IsoCon 0.3.3** (`pixi run -e isocon IsoCon pipeline`, defaults, `--nr_cores 4`) given all reads of the family (A + B + E,
FASTA, as Amendment 47's targeted input): its best final candidate against E. This is the baseline for the read-correction step, which
gets its own pre-registration.

### Reading, fixed now

- **S1, S2, S3 pass** → the projection is the stage's structure output: `<out>.structure.tsv` and `<out>.contigs.gtf` are written by
  default (they add files and change none), `o3_export.py` uses them, and the claim becomes "structure-correct in X of Y simulated
  reference-absent copies". The next level (real reads, bed H, annotated erased copies) is pre-registered before it runs.
- **S1 or S3 fails** → the projection is reported as a diagnostic, not a structure claim; the failing scenarios are named. No rule
  constant is changed in response; a new rule needs a new development round on new dev seeds.
- **S2 fails for some scenarios only** → structure is claimed for the passing scenarios only; the others are listed as limits.

### P0 result (2026-10-14, before any held-out world was built) — PASS

- Rust `candidates::structure` + `candidate_copies` phase 4 (binary sha1 bf1f968f, working tree on 49e62be3); prototype md5 unchanged
  (511562b8…); scorer `bench/o3_structure/score.py` md5 0a4315ef…, comparator `isocon_cmp.py` md5 6d7c141e….
- **Parity: 148 / 148 flagged dev candidates give the prototype's boundary list exactly** (s11, s12, s13 × E0, E2;
  `score.py parity`). One implementation defect was found and fixed on the way, in the implementation: the first port named the
  re-aligned segments `<candidate>|<start>_<end>` and gave 147 / 148. minimap2 breaks equally good hits by a hash of the read name
  (its `--seed`), so a terminal stretch with two equal hits (dev s12 E2 `cand_f12_0`, 2 mismatches on either copy) landed on the other
  copy. The segment names are therefore part of the reference behaviour; the binary now names them `<start>_<end>` as the prototype
  does, in rounds that never repeat a name within one minimap2 call.
- The seven earlier products are byte-identical to the pre-change binary's on all six dev runs (`cmp`).
- `cargo test --release`: 1,295 passed; the 13 failures are the known ones (3 `pipeline_driver_dispatch`, 10 missing fixtures on
  this laptop); the 9 new `candidates::structure` unit tests pass.

## OUTCOME (2026-10-14) — S1, S2, S3 all PASS on the held-out worlds

`pixi run bash bench/o3_structure/heldout.sh ~/rustle_data/o3_structure/heldout` (seeds 21-25 `--strands`, E0 and E2; binary sha1
bf1f968f): 270 absent copies, 242 flagged candidates, all 242 scored (each union's best true transcript is its own E).

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| S1 | structure-correct (±5), all but DIV@0.08 | 210 / 222 = 0.946 | ≥ 0.90 | ✅ PASS |
| S2 | LOSS / GAIN / DUP / SHUFFLE / TRUNC5 / INTRON | 30/30, 30/30, 17/20, **16/20 = 0.800**, 20/20, 19/20 | ≥ 0.80 each | ✅ PASS (SHUFFLE exactly at the bar) |
| S3 | boundary recall / precision (±5), all but DIV@0.08 | 1,553 / 1,561 = 0.995 / 1,553 / 1,561 = 0.995 | ≥ 0.97 each | ✅ PASS |

**Reported beside.**
- Exact (±0) structure: 185 / 222 = 0.833; GAIN 14/30 and SHUFFLE 10/20 exact (their novel-segment edges sit within a few bp), every
  other scenario at or near its ±5 count. DIV@0.08: 18 / 20 structure-correct. Strand: + 110/115, − 100/107.
- **The 14 misses, named:** 10 carry exactly one extra boundary and miss none (DUP s21 E2 f17, s24 f16 ×2; INTRON s22 E0 f22; DIV@0.02
  s23 E2 f04, s25 E0 f04; LOWEXPR s25 f25 ×2; DIV@0.08 s22 f09 ×2). The 4 SHUFFLE misses are all family f19 (s23, s24 × E0, E2), the
  layout with the absent copy on the strand OPPOSITE its home paralog (`strand_of`: k % 4 == 3): they miss 1-3 boundaries around the
  moved exons. No miss is a wrong placement of a recovered boundary.
- **Detection** (E flagged): DIV@0.005 0/20 (linked as an allele, by design: delta = 0.958%), DIV@0.01 12/20, every other scenario 20/20
  or 30/30, LOWEXPR (8-16 reads) 30/30. NULL families: no flag.
- **Reconstruction** (identity x coverage of the candidate against E ≥ 0.999), all 270 absent copies: **candidate_copies 242/270 = 0.896,
  IsoCon 256/270 = 0.948**. The whole difference is detection at ≤ 1%: DIV@0.005 0/20 vs IsoCon 17/20, DIV@0.01 12/20 vs 20/20.
  Wherever the copy is flagged, our candidate is exact (242/242); IsoCon is ≥ 0.999 in 239 of those 242 cases. IsoCon reports
  transcripts and has no test of reference absence; the 0.5% copies are what the allele cutoff delta excludes by design.

**Reading (as fixed).** The projection is the stage's structure output: `<out>.structure.tsv` and `<out>.contigs.gtf` are written by
default (new files; the earlier products are byte-identical). The claim: *in simulated families, `candidate_copies` gives a reference-
absent copy its exon structure within ±5 bp in 210 of 222 flagged cases (95%), across SNP divergence up to 4%, exon loss, gain,
duplication, shuffle, 5' truncation and resized introns, on both strands.* The named limit: shuffled exons on a copy oriented
opposite its home paralog. Not yet done: `o3_export.py` exports the unmapped-rescue chain's candidates and does not yet read
`contigs.gtf`; the next level (real reads, bed H, the erased copies' annotated structure) gets its own pre-registration on the main
machine. In simulation, reconstruction is not the bottleneck where a copy is detected (242/242 exact), unlike real reads (14/40 exact
on bed H), so the read-correction step must be developed against real-read error, not famsim's.


---

## PREREG_o3_family_best_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_family_best_2026-10-15.md`)*

## O3: a read counts only for the family holding its best alignment (2026-10-15)

Pre-registration, written before any `--psv-family-best` run and before the fresh held-out inputs were built. Code under
test: branch `o3-family-best`, commit `7c9d2e8b` (`copy_assign --psv-family-best`, opt-in). User goal (2026-10-15): "keep
only reads whose best alignment is in the family, or cluster related families together or require the unexplained alleles
to be carried by most of the cluster's reads". This test is the first option.

### Question

O3 calls a PSV cluster a copy missing from the genome when nothing in the genome explains it better. The real-data test
(register 1294, `docs/PREREG_o3_realdata_2026-10-15.md`) left `copy` clusters in about 30% of families whose copies were
all in the genome. Its follow-up (register 1295) traced most of them to reads from other loci, placed on the family through
secondary alignments by `--no-as-tied-only`.

For each read, the aligner already gives part of O3's test: a read whose best alignment lies at another locus is better
explained there. `--psv-family-best` lets a read into family F's clustering only when F holds one of its best alignments:
- a primary or secondary record carrying the read's highest AS (ties count wherever the primary went; supplementary records
  are segments, not placements);
- with an aligned base inside a copy of F.

Does this:
1. remove most `copy` clusters from intact families, and
2. keep the detection of copies removed from the genome?

### Known before running (recorded now)

- **Register 1295** classified the `copy` clusters of the added-back arms of heldout2 (57) and dev (27).
- **`primary_local` tallies** were computed while designing the rule, from the existing assignments. `primary_local` is the
  primary-only form of the rule, without ties. It is the share of a cluster category's reads whose primary alignment
  overlaps a copy of the family:

  | substrate | category | reads with `primary_local` |
  |---|---|---|
  | heldout2 | other family's copy | 34 / 2,292 |
  | heldout2 | not in the genome | 9 / 2,276 |
  | heldout2 | catalog lacks the exons | 165 / 320 |
  | heldout2 | own-copy variant | 27 / 129 |
  | heldout2 | reads in all other clusters | 16,829 / 26,169 |
  | dev | GWFAM175, absent from the reference | 469 / 469 |
  | dev | minority segment | 219 / 930 |
  | dev | catalog lacks the exons | 11 / 16 |
  | dev | own-copy variant | 4 / 39 |
  | dev | not in the genome | 19 / 59 |
  | dev | reads in all other clusters | 13,112 / 56,545 |

  So **heldout2 and dev are development substrates for this rule**. Their old-rule results are also known: recovered 4 / 12
  and 15 / 25; families with a `copy` cluster 17 / 60 and 16 / 53.
- **Without the flag**, this branch's `copy_assign` reproduces dev's added-back arm byte for byte (9 output files). The
  binary is frozen at `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`.
- **Not seen:** any `--psv-family-best` output, and which families the fresh held-out rule picks.

### Method (fixed)

**Fresh held-out families** (`heldout3`):
`prepare.py heldout --min-copies 2 --min-expressed 2 --exclude-picked .../heldout2/picked.tsv --seed 1 --cap 60`.
- Families of the shipped KB3781 catalog with >= 2 copies.
- At least 2 copies with >= 20 primary reads on their exons. Register 1294's rule required every copy; that is now relaxed.
- Every locus hull <= 200 kb, and no copy overlaps a copy of A13.
- Not one of heldout2's 60 families. The unrun 14-family build is not excluded: no output of it was ever made.
- The last expressed copy by (chrom, start) is masked. Silent copies stay in the catalog as survivors.
- A family whose masked hull overlaps any other catalog copy is dropped.
- More than 60 families: a seeded sample of 60.

**Reads, masking and re-mapping** are as in register 1294:
- up to 500 reads per copy;
- `bedtools maskfasta`;
- `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, the masked genome's own index, `GGO.splice.mmi`
  for the added-back arm.

**Runs.** Four on heldout3: {removed, added-back} × {old rule, new rule}. The old rule is register 1294's command:
`copy_assign --bam BAM --fasta GENOME --regions regions.txt --families COPIES.tsv --copies-fa COPIES.fa --no-as-tied-only
--psv-clusters --region-threads 1 --skip-poa-diagnostic --out PREFIX`. The new rule adds `--psv-family-best`. Both use the
frozen binary. On heldout2 and dev, only the new rule is run; their old-rule outputs exist.

### Scoring (fixed now; `bench/o3_realdata/score_fb.py`)

- A `copy` cluster is an `unexplained` cluster with est_div > delta = 0.00958.
- **Placed**, **eligible** (>= 6 masked-copy reads placed on the family) and **recovered** are register 1294's rules:
  - one `unexplained` cluster holds >= 0.8 of the masked copy's PLACED reads, with purity >= 0.9;
  - the denominator is the placed reads, which the filter does not change, so a read the filter drops counts against
    recovery;
  - the eligible set is the same under both rules.
- `cleared` is register 1294's rule.

| id | bar (fresh held-out, heldout3) |
|---|---|
| S1 | under the new rule, <= 0.10 of the families have a `copy` cluster in the added-back arm |
| S1b | the new rule's number of such families is <= 0.5 x the old rule's (paired, same families) |
| S2 | eligible masked copies recovered under the new rule >= those recovered under the old rule - 1 |
| G | denominator guard: families with an `explained` cluster in the added-back arm, new rule >= 0.9 x old rule |

- S2 needs >= 5 eligible masked copies. With fewer it is reported as not measurable, and the reading uses S1, S1b and G.
- S1b with 0 old-rule families is not measurable.
- G is there because a filter that deleted every read would clear every `copy` cluster. G checks that the families' own
  reads still form clusters the catalog explains.

### Readings (fixed now)

1. **S1, S1b, S2 and G pass.** Claimed, quoted with the counts: "keeping each read in the family that holds its best
   alignment removes most `copy` clusters of intact families on real gorilla reads, without losing detections". Recommended
   as the default inside `--psv-clusters`. That is the user's decision; the flag stays opt-in until then.
2. **S1, S1b and G pass, S2 fails.** Specificity gained at a cost in detection. Reported and not recommended; the lost
   detections are described.
3. **S1 or S1b fails, or G fails.** Not claimed. The remaining `copy` clusters are classified with register 1295's tools.

**Reported, not barred:**
- R1 and R2 of register 1294 under both rules;
- hosted (family, read) pairs kept by the filter;
- the masked copy's placed reads kept;
- `copy` clusters per family;
- the development substrates' numbers (heldout2, dev).

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-15)

Run with the frozen binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`, built from the Rust source committed as
`7c9d2e8b` (sha256 `b0eced3a99d93d4f...`), and `bench/o3_realdata/` on this branch. Outputs:
`/mnt/linuxdisk/tmp/o3_realdata/heldout3/`, plus `*_fb` and `score_fb.json` in `heldout2/` and `dev/`.

**Build.**
- `prepare.py` picked 40 families out of 326 candidates, after excluding heldout2's 60: 40 masked copies and 153 survivors.
  Reads: 7,510 masked-copy (D) and 15,692 survivor (S).
- Masked arm: 22,909 primary, 293 unmapped. Added-back arm: 23,202 primary.
- Cost: the masked index took 251 s at 22.6 GB; mapping, 8 parts of about 330 s; each `copy_assign` run, 226-293 s at
  <= 2.7 GB.

**The filter is the rule as specified.**
- `family_best.py`, an independent Python implementation, keeps exactly the pairs the Rust runs keep:
  - heldout3: 23,180 / 25,537 (added back) and 15,069 / 18,349 (removed);
  - heldout2: 26,793 / 31,798 (added back);
  - dev: 58,352 / 59,105 (added back);
  - the removed arms of heldout2 (8,649 / 12,135) and dev (46,347 / 49,090), so all six arms match.
- No clustered read is outside its kept set.
- The prototype re-clusters every remaining `copy` cluster at Jaccard 1.0.

| heldout3 (40 fresh families) | old rule | `--psv-family-best` | bar | |
|---|---|---|---|---|
| families with a `copy` cluster (added back) | 10 | 6 = 0.15 | S1 <= 0.10 | **FAIL** |
| the same, paired | 10 | 6 | S1b <= 5 | **FAIL** |
| `copy` clusters | 35 | 11 | | |
| families with an `explained` cluster | 36 | 36 | G >= 32.4 | pass |
| eligible masked copies recovered (24 eligible) | 7 | 5 | S2 >= 6 | **FAIL** |
| cleared once added back | 5 of 7 | 4 of 5 | | |

**Lost detections.** MCL117 (recall 0.879 -> 0.515) and MCL31 (0.841 -> 0.305). Once their copy is masked, their reads'
best alignment lies at a locus outside the catalog:
- every dropped masked-copy read of MCL117 (171) and MCL31 (48);
- 94 of heldout2 MCL126's 97 (the other 3 in MCL43 / MCL296).
So it is not a related family that holds them.

**Reported items (pre-registered).**
- The masked copy's placed reads kept by the filter: 2,428 / 3,391 (heldout3), 840 / 1,107 (heldout2), 5,252 / 7,286
  (dev).
- `copy` clusters per family, heldout3 added back, old -> new:
  - gone: MCL37 10 -> 0, MCL51 2 -> 0, MCL53 2 -> 0, MCL95 1 -> 0;
  - reduced: MCL32 9 -> 1, MCL63 2 -> 1;
  - unchanged: MCL45 5, MCL79 2, MCL125 1, MCL31 1.

| development | heldout2 (60) old -> new | dev (53) old -> new |
|---|---|---|
| families with a `copy` cluster | 17 -> 10 | 16 -> 16 |
| `copy` clusters | 57 -> 18 | 27 -> 27 |
| recovered | 4 -> 3 of 12 (MCL126 lost) | 15 -> 8 of 25 (7 lost) |

#### Verdict (2026-10-15)

**Reading 3: not claimed.**
- S1 and S1b fail: 6 of 40 families keep a `copy` cluster.
- S2 fails: 2 of 7 detections are lost.
- G passes.

On heldout3 the filter removes every cluster made of another family's reads (0 of the 11 left). But it costs detections: when a
copy is missing, its reads' best alignment often lies at a locus outside the catalog, an uncatalogued paralog.

The filter has a design limit. On heldout2, 3 new-rule clusters (MCL20 and MCL33, 2,894 reads) are still mostly other
families' reads. MCL20, MCL33 and MCL134 all have copies on the mitochondrion (`NC_011120.1`), so one best alignment there
counts as being in all three families.

**Known before running, corrected.** `primary_local` is always 0 for families assigned through `copy_assign`'s
cross-chromosome pass:
- dev, multi-chromosome families: 44,367 of 44,367 rows read 0, though the filter keeps them;
- one-chromosome families: it agrees with the filter on 14,076 of 14,091 rows.

So the tallies above understated what the rule keeps on dev. The column's doc comment says every per-family read rate must be
taken over it; that is a defect of its own, not fixed here.

#### The remaining `copy` clusters (register 1295's tools, as reading 3 requires)

heldout3 under the new rule: 11 clusters in 6 families, 270 reads. None is another family's reads.

| what the cluster is | clusters (reads) |
|---|---|
| chimpanzee transcripts, MCL45 (identity to chimp 0.997-0.998) | 2 (51) |
| MCL45, closest to chimpanzee (0.97) or nothing (0.85); MCL125 (0.92 to every ape) | 4 (89) |
| gorilla sequence outside both assemblies, MCL79 (0.994) | 1 (4) |
| the catalog lacks the exons: MCL63, MCL79 | 2 (66) |
| a minority segment, MCL31 | 1 (54) |
| an own-copy variant, MCL32 | 1 (6) |

#### Post hoc: would a variation graph with every isoform help? (user question; descriptive)

**More reference content, same exact rule** (`genomic_explain.py`). Register 1294's explained rule was asked against
three references for heldout3's 136 unexplained clusters (added back, old rule):

| reference | clusters explained |
|---|---|
| G0, each copy's spliced representative (as the method does) | 0 |
| G1, each copy's genomic locus +- 2 kb: every isoform, intron and PSV of the family | 1 |
| G2, G1 plus every locus where the cluster's reads align | 4 |

G1's one explanation is MCL0's own copy 3, an isoform its representative lacks, in an `allele_or_copy` cluster (est_div
0.0005). G2 adds three clusters explained by another family's copy: MCL51 twice by MCL107, and MCL95 by MCL23.

**What this proxy can and cannot say** (corrected after the branch review).
- G0, G1 and G2 all test one object: the cluster's consensus, built from its reads realigned to the HOST window. For reads
  from another locus, that consensus can sit far from where the reads really align. MCL31's cluster 0/1: the reads'
  median `de` at their own primary is 0.0007, yet the consensus is 0.088 from the host (est_div 0.053).
- In 28 of the 35 `copy` clusters the reads' median `de` is <= delta. Only 4 of those 28 have any G2 locus within delta of
  the consensus; the median distance is 0.13.
- So the proxy cannot separate the exact rule's strictness (one locus must carry all 20-338 alt sites) from the host-window
  realignment. Either can keep a cluster unexplained.
- What it does show: adding isoform content to this test explains almost nothing.
- The tolerant rule below works because it skips the host-window consensus and uses each read's own genome alignment.

**A tolerant rule** (`tolerant.py`). An unexplained cluster counts only when its reads' median best-alignment divergence
(minimap2 `de` at the primary, the aligner's best-scoring alignment anywhere in the genome; not necessarily the record with
the lowest `de`) exceeds delta:

| old-rule clusters, as run -> tolerant | heldout3 (40) | heldout2 (60) | dev (53) |
|---|---|---|---|
| families with a `copy` cluster | 10 -> **2** | 17 -> 6 | 16 -> 9 |
| `copy` clusters | 35 -> 7 | 57 -> 14 | 27 -> 17 |
| masked copies recovered | 7 -> 5 of 24 | 4 -> 3 of 12 | 15 -> 11 of 25 |

- Held out, the 2 remaining families are MCL45 (chimpanzee transcripts) and MCL79.
- Dev keeps GWFAM175, the genuine reference-absent copy.
- The lost recoveries are mostly `allele_or_copy` copies, already within delta of another locus (MCL22, MCL97, 3 of dev's
  4). MCL10's masked copy is labelled `copy`, but its reads still align at 0.9985 somewhere else in the masked genome:
  by O3's definition it is not absent.
- This was measured on substrates already seen, so it is a hypothesis for a new prereg on a fresh substrate. The gorilla
  pool is exhausted under these rules; human A119b is next.

**For a variation graph**, this says:
- The graph is the natural structure for the tolerant rule: a cluster is a missing copy when its consensus is more than delta
  from every path.
- In this proxy its content gain alone is small (G1 explains 1 of 136, with the caveat above). The gain comes from judging
  reads by their own best alignment, which the aligner already gives genome-wide.
- What a graph would add: both haplotypes, so allelic reads are explained; related families in one structure; and a clean
  path definition for the thesis.

#### Known minor issues (from the branch review; not fixed, no effect on these data)

- With no AS on any record, Rust's `BestAlignments` falls back to the LAST primary record, Python's `family_best.py` to
  the first.
- Rust scores a mapped record with an empty CIGAR, which Python skips.
- minimap2 output triggers neither.
- `score_fb.py` counts `copy_clusters` over every family in the run, not only the masked families. They are the same
  families on all three substrates.
- `genomic_explain.py`'s `g2_families` column names the NEAREST locus's families, which can differ from the locus that
  explains the cluster: MCL51's cluster 0/0 is explained by an MCL107 copy, while its nearest locus is outside the
  catalog.


---

## PREREG_o3_graph_real_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_real_2026-10-15.md`)*

## Graph-based O2 and O3 on handpicked real families (simulation v2, 2026-10-15)

Pre-registration, written before any held-out family's world was built. Code: branch `o3-graph`, commit `3dcf98a3`
(`bench/o3_graph/`).

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

This follows `docs/PREREG_o3_graph_sim_2026-10-15.md` (register 1297), whose worlds were famsim copies of single curated
genes in small genomes.

### Simulation v2 (`families.py`, `realworld.py`)

**Families.** Handpicked mainstay multi-copy families, taken by CAT gene-name pattern from CHM13 v2.0's CAT/Liftoff
annotation (the project's human default):
- copies linked when their canonical transcripts align at >= 90% identity (largest linked component);
- span <= 60 kb, non-overlapping, at most 10 copies (seeded);
- up to 4 curated isoforms per copy with distinct intron chains.

The split was fixed before any run:

| set | families (copies / isoforms) |
|---|---|
| development | SMN 2 / 8, NCF1 3 / 8, PAGE 3 / 8, CT45A 9 / 16, GAGE 10 / 12, POTE 10 / 28, GSTA 5 / 11, TP53TG3 4 / 13 |
| **held out** | NPIPA 6 / 17, NPIPB 10 / 24, TBC1D3 10 / 20, GOLGA8 10 / 24, GOLGA6L 4 / 5, CT47A 7 / 7, XAGE 3 / 12, SPATA31 7 / 14 |

NPIP, TBC1D3 and GOLGA8 are Soto's mainstays.

**Worlds.**
- The real CHM13 genome; reads mapped genome-wide with the baseline's minimap2 command.
- 60 reads per expressed copy over its isoforms: the canonical "main" isoform weighs 2.0, the others 0.5.
- HiFi errors from `bench/sim.py` (0.1% substitutions, 0.03% indels) and 0-30 bp end jitter.
- **Read names are a seeded shuffle and carry no truth** (review M5 of register 1297).
- The catalog's representative is each copy's main isoform.
- One SPECIAL copy per family (seeded) plays the scenario's part.

| scenario | what the special copy is |
|---|---|
| INTACT | an ordinary copy: every copy is in the genome and the catalog |
| REMOVED | hard-masked in the genome and dropped from the catalog; its reads still exist |
| OUTSIDE | in the genome, not in the catalog (an uncatalogued paralog) |
| BRANCH | the source of an extra ABSENT copy: its isoforms with a novel 150-bp exon at the middle of the main isoform's middle intron (spliced into every isoform with an intron across that point) and 0.2% exonic SNPs |

### Methods under test

**O2** (`o2graph.py`; the reads O2 places on the family = copy_assign's rows):

| method | what it does |
|---|---|
| current O2 | copy_assign `--no-as-tied-only`: rows with status `assigned` |
| **path O2** | each read aligned (minimap2 -x map-hifi, all secondaries) to EVERY curated isoform path of every catalogued copy. Distance to a copy = edits + unaligned bases against its best-fitting isoform. Assigned to the closest copy when the next is >= MARGIN farther; abstains when tied or thin, and as `unexplained` when even the closest copy is > delta away. **Primary: MARGIN = 2**; MARGIN = 1 reported. |
| path O2, representatives only | the ablation: each copy's representative path only (MARGIN 2) |

**O3** (`detect.py`; arms as in register 1297, plus S4):

| arm | what it does |
|---|---|
| CUR | `copy_assign --psv-clusters` |
| S123 | all-isoform explanation + branch calls + genome check |
| **S4** | the S1 clustering run only on the reads path O2 calls `unexplained`, plus the branch calls, kept by the genome check |

### Known before running (recorded now)

**Development results**, final code. O2, pooled over the four worlds (10,106 reads of catalogued copies, 1,449 of
non-catalogued ones):

| | assigned | correct | wrong | non-catalogued reads assigned |
|---|---|---|---|---|
| current O2 | 6,484 | 6,484 | 0 | 660 |
| path O2, margin 2 | 7,744 | 7,744 | 0 | 217 |
| path O2, margin 1 | 9,114 | 9,114 | 0 | 226 |
| path O2, representatives only, margin 2 | 6,502 | 6,348 | 154 | 170 |

O3 (8 families per scenario):

| | CUR | S123 | S4 |
|---|---|---|---|
| BRANCH detected | 0 | 8 | 8 |
| REMOVED detected | 3 | 3 | 4 |
| OUTSIDE false-call families | 4 | 0 | 0 |
| INTACT false-call families | 0 | 0 | 0 |

REMOVED's 8 removed copies include 4 that are detectable by O3's definition (the reads' median primary `de` on the masked
genome > delta: CT45A, GAGE, GSTA, PAGE). The 4 within delta are NCF1, POTE, SMN and TP53TG3; no arm calls them.

**Found and fixed on development:**
- copy_assign needs a family's span on a chromosome in one region, so regions are family hulls.
- Node-based graph O2 (allele nodes of a free GraphAligner alignment) makes artifactual calls. In GAGE it made 74 wrong of 289;
  the reads match both copies equally well linearly. Normalising the graph made it worse (173 wrong). Path-constrained
  alignment replaces it.
- A novel exon goes into every isoform whose intron spans it (not only the exact junction).
- The `novel` truth flag must look at both strands.
- S4 was added after GAGE's removed copy formed an impure cluster (59 of its reads with 48 of the host's, purity 0.55).
- BRANCH is scored on the absent copy's reads that carry the novel exon. Its other reads are the source copy's sequence
  within 0.2% SNPs, so no method can tell them apart; the all-reads result is reported beside it.

### Bars (held out: the 8 held-out families, seed 21; pooled over INTACT, REMOVED, OUTSIDE, BRANCH)

| id | bar |
|---|---|
| O2a | path O2 (margin 2) precision over reads of catalogued copies >= 0.995 |
| O2b | path O2 (margin 2) correct assignments >= 1.10 x the current O2's |
| O2c | path O2 (margin 2) assigns no more reads of non-catalogued copies than the current O2 |
| O2d | the representatives-only ablation makes more wrong assignments than all-isoform path O2 (step 1 is what keeps it precise) |
| O3a | S4 detects >= 6 of the 8 BRANCH absent copies |
| O3b | S4 detects every detectable REMOVED copy that CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| O3c | S4 has a false call in <= 1 of the 16 INTACT and OUTSIDE families |
| O3d | S4 detects >= 5 more BRANCH copies than CUR, and has no more false-call families than CUR |

Detectable (O3b) means the removed copy's reads have a median primary `de` on the masked genome > delta. That is O3's
definition: nothing in the genome explains the reads within delta.

### Readings (fixed now)

1. **O2a-O2d pass:** claimed, "on simulated reads of handpicked real families, O2 on every curated isoform path assigns
   >= 10% more reads correctly than the current O2 at >= 99.5% precision, and the representatives-only version is less
   precise".
2. **O3a-O3d pass:** claimed, "O2-informed O3 (S4) detects absent copies carrying a novel exon, keeps the current method's
   detections of sequence-divergent absent copies, and makes at most one false call over intact and uncatalogued-paralog
   families".
3. Each claim stands or falls separately. A failed bar is reported by name.

Scope either way: simulated reads (curated isoforms, HiFi errors) on the real genome. Real data's other species, fragments
and allele background are not modelled.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-15)

Run: code `549b5d1c`, the frozen copy_assign binary. Held-out worlds: the 8 held-out families, seed 21. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout/` (one directory per scenario, `o3_score.json`, `o2_*.json`).

**Disclosure: a scorer bug fixed after the first scoring.** Every held-out world carries the truth column `novel`, so
`score.py` took each REMOVED copy's source to be its novel-exon reads, which are none. That turned REMOVED's detections into
"no source" and its calls into false calls. The prereg scores BRANCH on novel-exon reads and REMOVED on all of the copy's
reads (register 1294), so the implementation was corrected to that rule (test added). Development's REMOVED world predates
the column and was unaffected. Only REMOVED's numbers changed.

**O2** (pooled over the four worlds: 12,766 reads of catalogued copies, 1,443 of non-catalogued ones):

| | assigned | correct | wrong | non-catalogued reads assigned |
|---|---|---|---|---|
| current O2 | 8,720 | 8,720 | 0 | 426 |
| **path O2, margin 2** | 9,964 | 9,964 | 0 | 181 |
| path O2, margin 1 | 10,650 | 10,650 | 0 | 213 |
| path O2, representatives only, margin 2 | 8,148 | 8,136 | 12 | 234 |

Per family (correct, current -> path O2 margin 2):

| family | reads | current | path O2 |
|---|---|---|---|
| GOLGA6L | 840 | 408 | 840 |
| NPIPA | 1,324 | 1,040 | 1,300 |
| NPIPB | 2,292 | 1,548 | 1,754 |
| TBC1D3 | 2,288 | 2,072 | 2,164 |
| GOLGA8 | 2,288 | 2,120 | 2,152 |
| SPATA31 | 1,564 | 1,532 | 1,564 |
| XAGE | 610 | 0 | 190 |
| CT47A | 1,560 | 0 | 0 (copies identical over their reads) |

**O3** (8 families per scenario):

| | CUR | S12 | S123 | **S4** |
|---|---|---|---|---|
| BRANCH detected | 0 | 8 | 8 | **8** |
| BRANCH false-call families | 3 | 3 | 0 | **0** |
| INTACT false-call families | 3 | 3 | 0 | **0** |
| OUTSIDE false-call families | 4 | 4 | 0 | **0** |
| REMOVED detected (2 detectable) | 0 | 1 | 1 | **1** |
| REMOVED false-call families | 3 | 3 | 0 | **0** |

- BRANCH holds on all of the absent copy's reads too: S4 8 / 8, CUR 0.
- 6 of the 8 removed copies are within delta of the masked genome by their reads' median `de`: CT47A 0.0016, NPIPB 0.0014,
  XAGE 0.0021, NPIPA 0.0026, SPATA31 0.0057, GOLGA6L 0.0093.
- The 2 detectable ones are TBC1D3 (0.0141; S4 detects it, CUR does not) and GOLGA8 (0.0143; missed by every arm).

| bar | value | |
|---|---|---|
| O2a | precision 1.0 | pass |
| O2b | 9,964 vs 8,720: +14.3% | pass |
| O2c | 181 vs 426 | pass |
| O2d | representatives only: 12 wrong; all isoforms: 0 | pass |
| O3a | BRANCH 8 / 8 | pass |
| O3b | 2 detectable removed copies (needs >= 3) | **not measurable** |
| O3c | 0 of 16 | pass |
| O3d | BRANCH +8; false-call families 0 vs 13 | pass |

#### Verdict (2026-10-15)

*Superseded by "Verdict, revised after the review" at the end of this document.*

**O2: reading 1, claimed** (simulated reads of handpicked real families). O2 on every curated isoform path assigns 14.3% more
reads correctly than the current O2, at precision 1.0. It also places fewer reads of non-catalogued copies (181 vs 426).
The representatives-only version makes wrong assignments (12 held out, 154 in development), so the all-isoform content is
what keeps it precise. Margin 1 assigns more still (10,650, no error) but is not the pre-registered primary.

**O3: claimed in part.** O3a, O3c and O3d pass. O2-informed O3 (S4):
- detects every absent copy carrying a novel exon (8 / 8; the current method 0);
- makes no false call over intact, uncatalogued-paralog, removed and branch families (0 / 32; the current method 13 / 32).

The REMOVED clause (O3b) cannot be tested on this set: 6 of its 8 removed copies are within delta of a remaining paralog, as
these young families' copies are. On the two that are not, S4 detects TBC1D3 and the current method neither.

**What the miss shows (post hoc).** GOLGA8's removed copy is a MOSAIC: its 58 reads form a pure cluster with 23 alt sites,
22 of them carried by some other copy, so no copy explains the cluster and too few sites are unexplained (`unresolved`). Its
reads are 1.43% from anything in the genome. A rule on the cluster's distance to its nearest single copy, not on
copy-private alleles, would call it. That is a method change for a new prereg on another family set.

**Scope:**
- simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome;
- real data's other species, fragments and allele background are not modelled;
- the special copies are one per family (seeded).

### Amendment 1 (2026-10-15): a catalog that is not copy-specific (the CROSS world)

Written after the whole-branch review and before the held-out CROSS world was built. Code: `9b6b98df`.

**Why.** The review found that much of path O2's gain is decided by the annotation rather than the sequence (review C1;
the corrections are below, under "Review corrections"). Every read here is simulated from one of its copy's own curated
isoforms, and each copy's catalog entry lists exactly those. Real annotation is not that copy-specific: an isoform
annotated on one paralog can be expressed by another. Then path O2 meets a read whose exon structure only a paralog's
catalog entry carries.

**World CROSS** (`realworld.py --scenario CROSS`):
- as INTACT: every copy is in the genome and the catalog, and its own-isoform reads are INTACT's;
- in addition, every copy expresses up to 2 isoform structures (seeded) annotated only on a paralog;
- such a structure is the paralog isoform's sequence spliced-aligned to the copy's own locus (minimap2 splice:hq). It is
  kept when the primary alignment covers >= 95% of the isoform with the same exon count, canonical junctions and an intron
  chain none of the copy's own isoforms has (junctions within 10 bp count as the same);
- the transcript is the copy's own genomic sequence over those exons, at the weight of a non-main isoform;
- the catalog still lists each copy's own isoforms only. Truth isoform `proj:<paralog>|<isoform>`.

**Known before running (development, seed 11).** 3,539 reads of catalogued copies, 771 of them from 67 paralog-annotated
structures over 46 copies:

| | assigned | correct | precision | paralog-structure reads: assigned / correct |
|---|---|---|---|---|
| current O2 | 2,228 | 2,219 | 0.996 | 400 / 391 |
| path O2, margin 2 | 2,459 | 2,223 | 0.904 | 303 / 67 |
| path O2, margin 1 | 2,864 | 2,620 | 0.915 | 361 / 117 |
| path O2, representatives only | 1,945 | 1,817 | 0.934 | 152 / 67 |

229 of path O2's 236 wrong calls name the paralog that carries the annotation. O3: every arm 0 / 8 false-call families.

**Bars (held out: the 8 held-out families, seed 21, world CROSS; wrong-family assignments count as wrong):**

| id | bar |
|---|---|
| X1 | path O2 (margin 2) precision over the reads of paralog-annotated structures >= 0.99 |
| X2 | path O2 (margin 2) precision over all reads of catalogued copies >= 0.995 (O2a's bar) |
| X3 | path O2 (margin 2) correct assignments >= the current O2's |
| X4 | S4 has a false call in <= 1 of the 8 families |

**Readings (fixed now):**
1. X1-X3 pass: O2's claim extends to a catalog that is not copy-specific.
2. X1 or X2 fails: O2's claim is scoped to a complete, copy-specific isoform catalog, and path O2 is not safe on real
   annotation as it is. A fix (every family structure threaded through every copy) needs its own prereg on a new family set.
3. X4 is read on its own: it pass/fails O3's specificity under such a catalog.

### Review corrections (whole-branch review, 2026-10-15; post hoc)

A fresh reviewer re-ran every number above and all of them reproduce. The prereg timing, the absence of leakage (anonymised
names; no method reads the truth) and the simulation's construction check out. The corrections below are post hoc; the bars
and their values above stand as computed.

**C1. About 42% of path O2's gain is decided by the annotation, not the sequence.**
- `o2graph.py seqsupport` aligns every hosted read to each catalogued copy's GENOMIC locus (spliced, +-1 kb). A read is
  genome-separated when its own copy's locus explains it at least MARGIN (2) edits better than any other copy's locus.
- 9,330 of the 12,766 reads of catalogued copies are genome-separated. On them path O2 makes 9,326 correct calls and the
  current O2 8,602: **+8.4%, below O2b's 10% bar.**
- Of path O2's 1,244 additional correct calls, 724 are genome-separated. The other 520 are reads another copy's locus
  explains within 1 edit; only the curated isoforms listed for each copy separate them. Every read is simulated from its
  copy's own curated isoforms, so those calls are correct by construction.
- Margin 1's extra calls are all of that kind (genome-separated correct: 9,326 at both margins).
- A read that aligns to one copy's isoform paths only is assigned without the margin test (`decide`); TBC1D3's gain is of
  this kind (reviewer).
- Representatives only (an incomplete catalog): 7,864 genome-separated reads correct, fewer than the current O2's 8,602.
- Hence Amendment 1 (above) and its results (below).

**I1. S4's O2 step adds nothing on the held-out set.** S4 equals S23 (branch calls + genome check, without the isoform step
or path O2) in every held-out family and scenario. Every BRANCH detection and the TBC1D3 REMOVED detection are branch calls;
every removed false call is removed by the genome check. Development's one difference (REMOVED GAGE, S4 only) did not
replicate. The O3 result belongs to the branch calls and the genome check, not to "O2-informed O3".

**I2. The current method's 13 false-call families are mostly another family's reads.**
- 12 of the 13 are GOLGA6L, NPIPA and NPIPB in each of the four worlds: reads of a sister family placed on them. Splitting
  families by gene-name pattern creates those sisters (NPIPA / NPIPB, GOLGA8 / GOLGA6L).
- Only TBC1D3 in OUTSIDE is the designed failure. There the current method scans the uncatalogued paralog only because the
  region hull includes the special copy (M2, below).
- The current method plus the genome check (S3) is also 0 / 32.
- copy_assign's opt-in `--psv-family-best` (register 1296) targets exactly such reads; it was not run here.
- So O3c / O3d's false-call comparison shows what the genome check does. It is not an advantage of S4 over a
  genome-checked current method.

**I3. Fixed:** S4 keyed path O2's unexplained reads by read only. A read unexplained in a sister family it was hosted in
entered its own family's S4 pool (about 550 reads per world). It is now keyed by (read, family), with a test. Re-run on all
four held-out worlds, every call is identical.

**I4. O2c (181 vs 426) comes mostly from BRANCH.** Reads of non-catalogued copies assigned, current vs path O2: BRANCH 240 vs 9,
REMOVED 119 vs 86, OUTSIDE 67 vs 86 (path O2 places more), INTACT 0 vs 0. In BRANCH the 150-bp novel exon puts the reads
past delta.

**M1 (graded up: a reported number was wrong). Fixed:** wrong-family assignments were not scored. The current O2 makes 56
(NPIPA 40, NPIPB 16), so its "wrong 0" was wrong: its precision is 0.9936. Path O2 makes none (1.0). Counting reads of
non-catalogued copies too: path O2 0.982, current 0.948.

**M3 (graded up: the O3b reading rests on it).** GOLGA6L's removed copy has median `de` 0.00935, against delta 0.00958. Had it
been above delta, O3b would have been measurable (3 detectable) and failed (S4 1 of 3).

**M4 (graded up: the only REMOVED detection).** S4's TBC1D3 detection is one branch call holding 50 of 60 source reads (0.83
against the 0.8 bar). The current method made two pure calls, of 17 and 42 reads, split across host copies; neither
reaches the bar.

**M5. Corrected:** CT47A's copies are not identical over their reads. At margin 1 path O2 assigns 480 CT47A reads correctly; at
margin 2, none.

**M7. Development provenance.**
- Development INTACT, OUTSIDE and REMOVED were built (18:59-19:17) before `realworld.py`'s final version (19:34).
- Development INTACT's path-O2 output predates the final `o2graph.py` by 24 s.
- Rebuilt with the final code, the three worlds' reads, truth, catalog and regions are byte-identical. Path O2 reproduces
  development INTACT's output exactly.

**Deferred (minor, not fixed):**
- M2: regions are built from every copy, the special one included, so copy_assign scans the OUTSIDE paralog (see I2).
- M6: the representatives-only run overwrites path O2's `graph/o2paths/*.tx.fa` intermediates; outputs are unaffected.
- M8: `families.py`'s tie-breaker does nothing.

#### Amendment 1 results (held out, world CROSS, 2026-10-15)

Run: code `9b6b98df`, prereg amendment `99932167`. 4,372 reads of catalogued copies; 940 of them come from 79
paralog-annotated structures over 57 copies (475 paralog isoforms projected). CT47A, GOLGA6L and XAGE (9 reads) have
almost none.

| | assigned | correct | wrong | precision | paralog-structure reads: assigned / correct |
|---|---|---|---|---|---|
| current O2 | 2,881 | 2,810 | 71 (all wrong-family: NPIPA 58, NPIPB 13) | 0.975 | 538 / 538 |
| path O2, margin 2 | 2,901 | 2,743 | 158 | 0.946 | 273 / 115 |
| path O2, margin 1 | 3,167 | 2,916 | 251 | 0.921 | 367 / 116 |
| path O2, representatives only | 2,250 | 2,232 | 18 | 0.992 | 95 / 83 |

- Genome-separated reads (3,036): correct, current 2,763 and path O2 2,535.
- Path O2 calls 1,350 reads `unexplained`; S4 forms no call from them once the genome check runs.
- O3: S4, S123, S23 and S3 have 0 false-call families of 8. CUR, R0, S1, S12 and S4n have 3 (GOLGA6L, NPIPA, NPIPB, the
  sister-family reads of I2).

| bar | value | |
|---|---|---|
| X1 | 115 / 273 = 0.42 | **fail** |
| X2 | 2,743 / 2,901 = 0.946 | **fail** |
| X3 | 2,743 vs 2,810 | **fail** |
| X4 | 0 of 8 | pass |

**Reading 2.** When a copy expresses a structure annotated only on a paralog, path O2 names the annotated paralog. Its
claim is scoped to a complete, copy-specific isoform catalog, and it is not safe on real annotation as it is.

#### Verdict, revised after the review (2026-10-15)

**O2: claimed only for a complete, copy-specific isoform catalog; not for real annotation.**
- O2a-O2d pass as computed, but only on the simulation's premise that each copy's catalog entry lists exactly the isoforms
  it expresses.
- Counting only reads the genome itself separates, path O2's gain is +8.4%, below O2b's 10%.
- When copies also express paralog-annotated structures (Amendment 1), path O2 is wrong on 58% of the reads it assigns
  from them. Its overall precision drops to 0.946, and it assigns fewer reads correctly than the current O2 (X1-X3 fail).
- With representatives only (an incomplete catalog), it gets fewer genome-separated reads right than the current O2.
- The current O2's precision is 0.9936 with wrong-family calls counted (0.975 in CROSS; all of them sister-family calls).
- Next, under its own prereg on a new family set: thread every family structure through every copy (each copy's paths =
  its own isoforms + every paralog isoform projected onto its locus), so that only the sequence separates copies. Its
  ceiling on this set is the genome-separated gain (+8.4%).

**O3: the branch calls and the genome check, not O2-informed O3.**
- BRANCH 8 / 8 (current 0) holds, from the branch calls.
- No false call over 32 families (and 0 / 8 in CROSS), from the genome check. The current method plus the genome check
  does as well; the current method's own false calls are mostly sister-family reads.
- REMOVED is not measurable: 6 of 8 removed copies are within delta, one of them by a hair (0.00935 vs 0.00958).
  TBC1D3's detection is marginal (0.83 vs 0.8).
- The GOLGA8 mosaic rule remains the next O3 step (new prereg).

Scope unchanged: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome.


---

## PREREG_o3_graph_round4_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round4_2026-10-15.md`)*

## Graph O2 as assign-or-abstain, and family-pooled mosaic-aware O3, on a fifth handpicked family set (2026-10-15)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

### Where rounds 1-3 left it (registers 1298-1301, after their reviews)

No O2 or O3 claim has survived as pre-registered.
- **O2.** Graph O2's gains in correct assignments were decided by the annotation: per-copy isoform lists (round 1), one
  family's non-canonical intron (FCGR, round 2), and annotated end positions (round 2's development families). With
  padded paths, nothing was gained on round 3's fresh set (x0.995).
- **O2, what held every time:** graph O2 keeps the current O2's correct calls within 2% at precision >= 0.999. It also
  abstains on reads no catalogued copy explains within delta, placing far fewer reads of non-catalogued copies (round 2:
  223 vs 608; round 3: 375 vs 891). An annotation-free locus O2 abstains the same way, but it loses 5-6% of the correct
  calls.
- **O3.** Branch calls and the genome check were right every time: 8 / 8 novel-exon absent copies on three fresh sets, no
  false call after the genome check. Removed copies failed the bar twice.
  - Round 2: PSG's reads were split over four hosts.
  - Round 3: KIR was missed because the mosaic rule was coded more narrowly than documented. It skipped linked clusters
    with few copy-private alleles; fixed in `b0a3e59f`.

This round tests the two things that held or were fixed, on a fresh family set, as one claim each.

### Methods under test

**O2: padded union O2** (`o2graph.py run WORLD unionp 2`), unchanged from round 3. Each copy's paths are its own curated
isoforms plus every paralog isoform spliced-projected onto its locus, with both ends padded by 50 bp of the copy's genome.
A read is assigned to the closest copy at margin 2 and abstains (`unexplained`) beyond delta. Beside it: the current O2
(copy_assign `--no-as-tied-only`) and locus O2 (the same decision on each copy's genomic locus).

**O3: H** (`detect.py`), as round 3, with two fixes from the round-3 review (`b0a3e59f`):
- the mosaic rule, as documented: every cluster no copy explains is called when its nearest single copy is > delta away;
- one call per pooled member set.

H clusters the reads padded union O2 cannot explain, pooled over all of the family's hosts. The mosaic rule and the
branch calls apply, and the genome check keeps a call only when the genome does not explain its reads.

### Worlds

As round 3, seed 51: INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon), and CROSS
(`--projection liftover`). Worlds share INTACT's reads of catalogued copies (round-2 review M4), so the O2 bars are set per
world.

### Families (held out, picked before any run)

Candidates in a fixed order: TRIM, LILR, DEFB10, CSAG, MAGEB, CTAG, SERF1, CBWD, MAGEC, CYP2A, NBPF, GTF2I. Too few qualified,
so more candidates were added in a fixed order before any run: PRAMEF, VCX, SPDYE, FAM86, GPR89, STAG3L, POM121, CCDC144. The
set is the first 8 with >= 3 copies, so that a removed copy leaves >= 2 catalogued hosts (round 3's CYP4F had none).

| family | copies / isoforms |
|---|---|
| TRIM (TRIM49, 49B, 49C) | 3 / 5 |
| DEFB10 (DEFB103A/B) | 3 / 3 |
| CSAG | 4 / 13 |
| CYP2A (2A6, 2A7, 2A13) | 3 / 9 |
| NBPF (7, 8, 9, 11, 12) | 5 / 15 |
| GTF2I (pseudogenes GTF2IP2, P3, P8) | 3 / 3 |
| PRAMEF | 9 / 11 |
| VCX | 4 / 9 |

`/mnt/linuxdisk/tmp/o3_graph/families4.json` (earlier sets identical to `families3.json`; `families4_all.json` holds every
candidate).

### Known before running (development: the four earlier family sets)

Development sets: dev (the original development families), dev-2 (round 1's held-out set), dev-3 (round 2's), dev-4
(round 3's), all seen now.

O2, the bars on development (padded union O2 against the current O2 and locus O2):

| set | A1: lowest world ratio of correct calls | A2: lowest world precision | A3: non-catalogued placed, padded union vs current (OUTSIDE, REMOVED, BRANCH; pooled ratio) | A4: vs locus O2 |
|---|---|---|---|---|
| dev | 1.017 | 1.0 | 118 vs 118, 118 vs 228, 8 vs 314; 0.37 | x1.041 |
| dev-2 | 1.083 | 1.0 | 87 vs 67, 87 vs 119, 0 vs 240; 0.41 | x1.125 |
| dev-3 | 1.060 | 1.0 | 93 vs 119, 93 vs 162, 37 vs 327; 0.37 | x1.053 |
| dev-4 | 0.983 (CROSS) | 1.0 | 183 vs 241, 183 vs 373, 9 vs 277; 0.42 | x1.052 |

In OUTSIDE the uncatalogued paralog is in the genome and mostly within delta of a catalogued copy, so no abstention can
tell its reads apart (by O3's own definition they are explained). Padded union O2 places as many as the current O2 on dev
(118 vs 118) and more on dev-2 (87 vs 67). A3's "fewer" clause therefore covers REMOVED and BRANCH; OUTSIDE is reported.

O3, H with the post-review fixes (`b0a3e59f`), re-run on dev-3 and dev-4:

| set | BRANCH detected | REMOVED detected / detectable | false-call families, INTACT + OUTSIDE + CROSS | Hn false-call families (before the genome check) |
|---|---|---|---|---|
| dev-3 | 8 / 8 | 3 / 3 (FCGR, LRRC37A, PSG) | 0 / 24 | 8 world-families |
| dev-4 | 8 / 8 | 3 / 4 (CGB, KIR, UGT2B; CYP4F unreachable) | 0 / 24 | 6 world-families |

The same re-run on dev and dev-2 was still running when this was committed. Its numbers are added under Results; they do
not change any bar.

### Bars (held out: the 8 families above, seed 51; wrong-family calls count as wrong)

| id | bar |
|---|---|
| A1 | padded union O2 correct >= 0.98 x the current O2's, in every world |
| A2 | padded union O2 precision >= 0.995 over reads of catalogued copies, in every world |
| A3 | in REMOVED and in BRANCH padded union O2 places fewer reads of non-catalogued copies than the current O2, and <= 0.6 x as many over OUTSIDE, REMOVED and BRANCH together |
| A4 | padded union O2 correct >= 1.03 x locus O2's, pooled over the five worlds |
| B1 | H detects >= 6 of the 8 BRANCH absent copies |
| B2 | H has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| B3 | H detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| B4 | H detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |

Detectable (B3): the removed copy's reads have a median primary `de` on the masked genome > delta.

### Readings (fixed now)

1. **A1-A4 all pass:** claimed, "on simulated reads of handpicked real families, O2 on a family graph that threads every isoform
   structure through every copy keeps >= 98% of the current O2's correct assignments at >= 99.5% precision. It abstains on
   reads no catalogued copy explains, so it places fewer reads of non-catalogued copies (<= 60% as many). It keeps >= 3% more
   correct assignments than an annotation-free locus O2 with the same abstention." Any failing bar: not claimed.
2. **B1-B4 all pass:** claimed, "graph O3 that clusters a family's O2-unexplained reads together, with the mosaic rule, detects
   novel-exon absent copies and >= 75% of detectable removed copies, with at most one false call over intact,
   uncatalogued-paralog and paralog-annotation families". B3 not measurable (< 3 detectable): claimed without the
   removed-copy clause. B3 measurable and failing, or any other bar failing: not claimed.
3. Each claim stands or falls as a whole. A failed bar is reported by name.

Reported beside the bars:
- path O2;
- G, S23, S123m and Hn (H without the genome check);
- O2 per family;
- the liftover CROSS structures whose intron chain union does not thread (round-2 review I2).

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome. Real data's other species,
fragments and allele background are not modelled.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (held out, 2026-10-15)

Run: code and prereg `f788b838` (the detector as of `b0a3e59f`), the frozen copy_assign binary, seed 51. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout4/` (`o2r4_*.json`, `o3_score.json`).

**O2** (correct / assigned; there were no wrong-family calls):

| world | reads | current | padded union O2 | A1 ratio | locus O2 | path O2 | non-catalogued placed: current / union / locus |
|---|---|---|---|---|---|---|---|
| INTACT | 2,047 | 1,855 / 1,855 | 1,850 / 1,850 | 0.997 | 1,736 | 1,850 / 1,850 | - |
| CROSS (liftover) | 2,547 | 2,325 / 2,325 | 2,252 / 2,252 | **0.969** | 2,207 | 2,100 / 2,115 | - |
| OUTSIDE | 1,566 | 1,434 / 1,434 | 1,431 / 1,431 | 0.998 | 1,376 | 1,431 / 1,431 | 60 / 60 / 60 |
| REMOVED | 1,566 | 1,437 / 1,437 | 1,431 / 1,431 | 0.996 | 1,376 | 1,431 / 1,431 | 121 / 60 / 60 |
| BRANCH | 2,047 | 1,855 / 1,855 | 1,850 / 1,850 | 0.997 | 1,736 | 1,850 / 1,850 | 324 / 28 / 28 |
| pooled | | 8,906 | 8,814 (x0.990) | | 8,431 | | 505 / 148 / 148 |

CROSS's 500 paralog-structure reads:
- Padded union O2 gets 402 / 402 right and abstains on the rest: VCX 42, CSAG 34, NBPF 9 and CYP2A 1 `unexplained`, PRAMEF 12
  `thin`.
- The current O2 gets 470 / 470.
- Path O2 gets 250 / 265.
- The liftover structures that union does not thread cost it 68 reads the current O2 places correctly.

**O3** (8 families per world; H with the post-review detector):

| | CUR | S23 | S123m | **H** |
|---|---|---|---|---|
| BRANCH detected | 1 | 7 | 7 | **7** |
| INTACT + OUTSIDE + CROSS false-call families (of 24) | 5 | 0 | 0 | **0** |
| REMOVED detected (6 detectable) | 5 | 5 | 6 | **6** |
| REMOVED + BRANCH false-call families | 0 | 0 | 0 | **0** |

- **Detectable removed copies** (median `de` on the masked genome): NBPF 0.096, CYP2A 0.036, GTF2I 0.022, PRAMEF 0.021,
  TRIM 0.022 and VCX 0.015. CSAG (0.0039) and DEFB10 (0.0) are within delta.
- **PRAMEF** is detected only by the arms with the mosaic rule (S123m, H).
- **The BRANCH miss** (every graph arm) is CYP2A: 44 source reads, one call short of the recovery rule.
- **G** has no union O2 table in this round (round 4 ran only padded union O2), so its PSV calls are void here and it is not
  reported.
- **Hn**, before the genome check, has false calls in CROSS (NBPF, VCX) and OUTSIDE (CYP2A, GTF2I, NBPF, PRAMEF, TRIM, VCX).
  The genome check removes all of them.

| bar | value | |
|---|---|---|
| A1 | lowest world 0.969 (CROSS) | **fail** |
| A2 | 1.0 in every world | pass |
| A3 | REMOVED 60 < 121, BRANCH 28 < 324; pooled 148 / 505 = 0.29 | pass |
| A4 | 8,814 vs 8,431: x1.045 | pass |
| B1 | 7 / 8 | pass |
| B2 | 0 of 24 | pass |
| B3 | 6 of 6 detectable (1.0), including CUR's 5 | pass |
| B4 | H 13 vs S23 12, false-call families 0 = 0 | pass |

**Development re-runs finished after the commit** (H, post-review detector), as the prereg said they would be reported:
- dev: BRANCH 8 / 8, REMOVED 4 / 4, no false call;
- dev-2: BRANCH 8 / 8, REMOVED 2 / 2, but a false call in NPIPB in INTACT, CROSS and OUTSIDE (3 of 24) and in REMOVED. In each,
  the call is a pure cluster of NPIPA.5's reads placed on the sister family NPIPB. The mosaic rule calls it, since NPIPA.5 is
  no NPIPB copy. The genome check keeps it because those reads' median `de` is 0.0097, just above delta (0.00958).

So the mosaic rule's false-call guard is the genome check alone, and it fails when another family's reads land on a family
and sit at delta. Round 4's held-out set had no sister families.

#### Verdict (2026-10-15)

**O3: claimed (reading 2; B1-B4 pass).** "On simulated reads of handpicked real families, graph O3 that clusters a family's
O2-unexplained reads together, with the mosaic rule:
- detects novel-exon absent copies (7 of 8; the current method 1);
- detects >= 75% of detectable removed copies (6 of 6; the current method 5);
- makes no false call over intact, uncatalogued-paralog and paralog-annotation families (0 of 24; the current method 5)."

Scope and caveats:
- The removed-copy gain over the current method is the mosaic rule (PRAMEF). Here pooling detected nothing per-host
  clustering with the mosaic rule did not (S123m = H).
- Development set 2 shows the failure mode above: sister-family reads at delta.
- This is the first O3 claim of four rounds; rounds 2 and 3 failed their removed-copy bar with the code as then frozen.

**O2: not claimed (A1 fails in CROSS, 0.969 vs 0.98).** Padded union O2 kept >= 99.6% of the current O2's correct calls in the
four worlds where every expressed structure is threaded, at precision 1.0. It placed 71% fewer reads of non-catalogued copies
in REMOVED and BRANCH, and kept 4.5% more correct calls than locus O2. But where copies express paralog structures it does
not thread, it abstains on reads the current O2 places correctly (68 of 500 in CROSS). Across four fresh sets, graph O2
never gained correct assignments that held up: the gains were annotation-decided or absent.

Scope: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome.

### Review corrections (round-4 review, 2026-10-16; post hoc)

A fresh reviewer re-derived every round-4 number from the raw outputs, and all of them reproduce. Timing is clean: prereg
23:28:40, first held-out file 23:28:53. No Critical finding. The corrections are about attribution and scope.

**Credit (I1).**
- The false-call drop is the genome check's: S3, the current method plus the genome check with no graph, also has 0 of 24
  false-call families and detects 5 of 6 removed copies.
- The novel-exon detections are the branch calls'.
- The one removed copy H finds beyond the current method, PRAMEF, comes from the mosaic-rule extension in `b0a3e59f`. Its
  cluster is `unexplained` at est_div 0.0076-0.0089, below delta, and its nearest copy is 0.019-0.023 away.
- Pooling added nothing here: S123m equals H in every world.
- The pre-registered sentence says "at most one false call", not "no false call".

**Scope (I2).**
- On development set 2 the claimed detector fails B2: NPIPB false calls in INTACT, CROSS and OUTSIDE (3 of 24), plus REMOVED.
  `b0a3e59f` introduced them; H before the fix had none there.
- Development set 2 is the only one of the five sets with sister-family hosting: 550-790 reads per world are hosted on
  another family, against 0 in this held-out set and in dev-3 and dev-4. So this held-out 0 / 24 does not test that case.
- Part of it is a simulation artifact. NPIPA.5's Liftoff models put every junction 2-5 bp off a canonical site, so its reads
  sit at median `de` 0.0092 even on their own locus. The 10-11-read clusters reach 0.0097, above delta by 0.00012.

**Open code gaps, not fixed (I3).** Either fix is a method change that needs its own prereg; no round-4 outcome depends on
them.
- `status_of` returns `explained` before any distance is computed, and `explains` compares a copy only at sites called from the
  reads. A gene-conversion-type absent copy (a host plus a tract of another copy's alleles) comes out `explained` although
  every single copy is 1-2% away.
- `nearest_copy` counts deletions shared by the copy and the reads as differences, because sites are called from A/C/G/T only.
  KIR.6's own reads come out 0.0135 from their own copy. So the mosaic distance has a noise floor near delta.

**Reproducibility (I4).** B3's detectable set and A1-A4 came from scratch scripts. They are now `bench/o3_graph/summary.py`
(`o3`, `o2bars`; tests), which reproduces every value above.

**O2 history (I5).** The sentence "graph O2 never gained correct assignments that held up" overstated it. Round 1's
genome-separable gain was +8.4%, failing its 10% bar. Round 2's was about 1.5%, round 3's none, and round 4's -1%.

**Smaller corrections.**
- **m1.** The BRANCH miss (CYP2A) is not "one call short". One branch call holds all 44 scored reads plus 8 reads of the same
  absent copy without the whole novel exon, so purity is 0.846 against 0.9. `detected_all_reads` is true.
- **m2.** "71% fewer" is the three-world figure (148 / 505). REMOVED and BRANCH alone are 88 / 445, 80% fewer. Locus O2
  places exactly as few (60 / 60 / 28), so A3 measures the shared delta abstention, not the graph.
- **m3.** 15 of the 41 lifted CROSS structures are not threaded (158 reads, 73 still right). The 68 lost reads are 67 of
  those plus 1 from a threaded CYP2A structure.
- **m4.** Dropping OUTSIDE from A3's "fewer" clause before commit was disclosed and data-free, but the stated reason does not
  hold. Padded union calls most of the uncatalogued paralog's reads `unexplained` (dev 304 / 483, dev-2 218 / 481, held out
  363 / 481); both methods place only the minority within delta. The change decided A3 here (OUTSIDE 60 = 60 would have failed
  the original wording), but not the verdict.
- **m5.** "Precision >= 0.999 every time" is false for round 1 (CROSS 0.946) and for round 2 per world (CROSS 0.9973). The
  final detector's record across the earlier sets: dev-4 sits exactly at B3's 0.75 (CGB `de` 0.0101), and dev-2 fails B2.
- **m7.** The duplicated "## Results" heading is kept as in the earlier docs: the first marks where the prereg ends.

**Every H detection** (recall, purity, genome `de`):
- REMOVED: CYP2A 61/61, GTF2I, NBPF, PRAMEF, TRIM and VCX 60/60 (VCX also 58/60), purity 1.0; de 0.015-0.096.
- BRANCH: CSAG 40/40, DEFB10, GTF2I, NBPF, PRAMEF, TRIM and VCX 60/60.
- No detection sits at the 0.8 recall bar (the lowest is 0.967). Losing PRAMEF would still pass B3 and B4.

#### Verdict, revised after the review (2026-10-16)

**O3: claimed (reading 2; B1-B4 pass on the held-out set), in the pre-registered words.** "On simulated reads of handpicked
real families, graph O3 that clusters a family's O2-unexplained reads together, with the mosaic rule, detects novel-exon
absent copies and >= 75% of detectable removed copies, with at most one false call over intact, uncatalogued-paralog and
paralog-annotation families." Observed: 7 of 8, 6 of 6, 0 of 24.

Attribution and scope:
- the novel exons come from the branch calls;
- the false-call drop comes from the genome check (the current method plus it, S3, also gets 0 / 24 and 5 / 6);
- the extra removed copy (PRAMEF) comes from the mosaic-rule extension;
- pooling added nothing on this set;
- the claim holds for families without sister-family hosting. On development set 2, where another family's reads land on
  the family near delta, this detector makes false calls in 4 of 5 worlds.

**O2: not claimed (A1 fails, 0.969 in CROSS).**
- Padded union O2 is precise: 1.0 in every world here.
- It places as few reads of non-catalogued copies as locus O2, which shares its delta abstention.
- It keeps 4.5% more correct calls than locus O2.
- It abstains on paralog structures it does not thread.
- Across five family sets its genome-separable gain over the current O2 has been +8.4% at most (round 1, failing its bar),
  and about 0 since.


---

## PREREG_o3_graph_sim_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_sim_2026-10-15.md`)*

## O3 with a variation graph: three steps, tested on curated ideal scenarios (2026-10-15)

Pre-registration, written before any held-out world was built. Code: branch `o3-graph`, commit `883e54bc` (`bench/o3_graph/`).

**User goal (2026-10-15):** "proceed with your 3 step recommendation for giving the vg a bit more detection power; ... we
can use curated isoforms from refseq or cat to simulate ideal scenarios".

The recommendation followed the two-locus graph spike (`docs/PREREG_o3_family_best_2026-10-15.md`, Results, and memory).
It was: feed the graph into the current PSV-cluster method (`copy_assign --psv-clusters`, prototype
`bench/o3_unexplained/psvclust.py`) rather than build a separate detector.

### The three steps (as arms of `bench/o3_graph/detect.py`)

| arm | what it adds |
|---|---|
| CUR | the current method: `copy_assign --no-as-tied-only --psv-clusters`; a call = an `unexplained` cluster with est_div > delta |
| R0 | the prototype's per-host clustering against each copy's representative: CUR reproduced in Python |
| S1 | **step 1**: a cluster is explained against each copy's EVERY curated isoform, the isoform paths of the family's variation graph, not its one representative |
| S12 | S1 + **step 2**: branch calls. GraphAligner maps the family's hosted reads to its variation graph (`vg msga` of its copies' loci, every curated isoform added as a path with `vg rna`). A read with >= 50 bases on no graph path is a branch read: large insertions, dense error windows where the aligner forced novel sequence onto the graph, gaps, unaligned ends. Reads sharing a novel segment are linked; a group of >= 5 reads is a call. |
| S123 | S12 + **step 3**: a call is kept only when the genome does not explain its reads (see below) |
| S3 | R0 + step 3 alone |

Step 3 keeps a call only when the genome does not explain its reads:
- a PSV call: the median minimap2 `de` of its reads' primary alignments is > delta;
- a branch call: its reads' median unexplained bases on the genome (soft clips and dense error windows) are >= 50, or that
  median `de` is > delta.

delta = 0.00958, as in registers 1294-1296.

### Ideal scenarios (`bench/o3_graph/world.py`, famsim)

**Templates** (`templates.py`): curated CAT genes on CHM13 chr20 (`bakeoff/human_chr20/chr20_ref_cat.gtf`).
- Each has 4-14 exons, a span <= 40 kb, every intron GT..AG, and is not a transcript-assembly `MSTRG.*` gene.
- Each has >= 1 curated isoform that is an exon-skipping variant of the canonical and skips an internal cassette exon.
- 120 genes qualify. A seeded shuffle (seed 1) picks 24: the first 12 for development, the next 12 held out.
  - Development: ZNF217, SALL4, C20orf197, SLCO4A1-AS1, NSFL1C, GSS, AL031674.1, RBM38, FAM182A, AL160412.1, ARFGAP1, SNX5.
  - Held out: C20orf27, RBCK1, YTHDF1, POFUT1, PPDPF, LINC01370, SNRPB, AL139350.1, SAMD10, ZBP1, SPAG4, SNHG17.

**Copies and isoforms.**
- Every copy expresses the template's curated isoforms: the canonical, the MAIN isoform (skips the first cassette exon) and
  up to two more.
- The main isoform is the most expressed, so it is the catalog's representative, which lacks the cassette exon as real
  representatives do.
- 60 reads per expressed copy; famsim errors 0.1% substitution, 0.03% indel.

**Scenarios.** One world per scenario, one family per template; A = the template, B = A + 3% SNPs.

| scenario | copies | what should happen |
|---|---|---|
| NULL | A, B | no call |
| ALLELE | A, B, H = A + 0.4% (absent) | no call (an allele) |
| ABSENT@0.01 / @0.02 | A, B, E = A + 1% / 2% (absent) | E is called |
| BRANCH | A, B, E = A + a novel 150-bp exon + 0.2% SNPs (absent) | E is called: structure differs, sequence barely does |
| OUTSIDE@0.01 / @0.03 | A, B in the catalog, D = A + 1% / 3% in the genome only | no call: the genome explains D |
| TWIN_ALT | A, T = A + 15% SNPs in the cassette exon only, B | no call: T is a present copy |

**Scoring** (`score.py`, register 1294's rules):
- **Detected:** one call holds >= 0.8 of the absent copy's placed reads at purity >= 0.9.
- **False call:** a call whose reads are less than half the absent copy's. In a family with no absent copy, every call is
  false: NULL, ALLELE, OUTSIDE, TWIN_ALT.

### Known before running (recorded now)

**Development** (12 templates per scenario, seed 11), final code:

| scenario | CUR | R0 | S1 | S12 | S123 | S3 |
|---|---|---|---|---|---|---|
| ABSENT@0.01 detected | 8 | 8 | 8 | 8 | 7 | 7 |
| ABSENT@0.02 detected | 12 | 12 | 12 | 12 | 12 | 12 |
| BRANCH detected | 0 | 0 | 0 | 12 | 12 | 0 |
| OUTSIDE@0.01 / @0.03, families with a false call | 8 / 12 | 8 / 12 | 8 / 12 | 8 / 12 | 0 / 0 | 0 / 0 |
| TWIN_ALT, families with a false call | 4 | 4 | 0 | 0 | 0 | 0 |
| NULL, ALLELE, families with a false call | 0 | 0 | 0 | 0 | 0 | 0 |

The ABSENT@0.01 family S123 loses (SALL4) sits at delta: est_div 0.00979, the reads' median `de` 0.0095.

**Fixes made on development, before this prereg:**
- one world per scenario: in a mixed world, the same template's families are copies of each other;
- the novel exon goes into a long intron;
- unaligned ends, and dense error windows in GraphAligner's and minimap2's alignments, count as novel sequence. Both
  aligners force a random insert onto the reference as small mismatches and indels.
- novel segments are linked with short-sequence minimap2 parameters (the asm presets need an alignment score >= 200).

The current method and R0 make the same calls on development (the prototype = Rust, as on real data).

### Bars (held out: templates 13-24, seed 21; 12 families per scenario)

| id | bar |
|---|---|
| H1 | S123 detects >= 11 / 12 ABSENT@0.02 and >= 10 / 12 BRANCH, and ABSENT@0.01 >= CUR's count - 2 |
| H2 | S123 has a false call in <= 3 of the 60 families with no absent copy (NULL, ALLELE, OUTSIDE x2, TWIN_ALT) |
| H3 | S123 detects >= 8 more BRANCH families than CUR, and has >= 10 fewer false-call families than CUR over those 60 |
| H4 | each step is needed where it should be: (a) S1 has fewer TWIN_ALT false-call families than R0 (not measurable if R0 has none); (b) S12 detects >= 8 more BRANCH families than S1; (c) S123 has >= 10 fewer OUTSIDE false-call families (both rates pooled) than S12 |

### Readings (fixed now)

1. **H1-H4 pass.** Claimed, quoted with the counts: "In ideal scenarios built from curated isoforms, the graph-informed O3
   detector finds absent copies that differ only in structure, which the current method cannot. It also stops false calls
   from uncatalogued paralogs and from copies a representative makes look identical, and keeps the current detection of
   sequence-divergent absent copies; each step is needed for its case." Scope: simulation; real data not yet tested.
2. **H1 and H2 pass, H3 or H4 fails.** The detector works; the failing comparison is reported as such.
3. **H1 or H2 fails.** Not claimed.

Reported, not barred: every arm's counts per scenario; recall and purity of the calls; S3 (step 3 alone).

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-15)

Run: code `23441fda` (as `883e54bc`; the prereg adds no code). Held-out worlds: templates 13-24, seed 21. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/heldout/` (one directory per scenario, `score.json`). Each world takes about 30 s: build, map,
the current method, detector.

| held out (12 families each) | CUR | R0 | S1 | S12 | **S123** | S3 |
|---|---|---|---|---|---|---|
| ABSENT@0.01 detected | 8 | 8 | 8 | 8 | **8** | 8 |
| ABSENT@0.02 detected | 12 | 12 | 12 | 12 | **12** | 12 |
| BRANCH detected | 0 | 0 | 0 | 12 | **12** | 0 |
| OUTSIDE@0.01 / @0.03, families with a false call | 8 / 12 | 8 / 12 | 8 / 12 | 8 / 12 | **0 / 0** | 0 / 0 |
| TWIN_ALT, families with a false call | 9 | 9 | 0 | 0 | **0** | 0 |
| NULL, ALLELE, families with a false call | 0 | 0 | 0 | 0 | **0** | 0 |

- Families with a false call, of the 60 with no absent copy: CUR 29, S123 0.
- S123's detected calls: median recall 1.0, purity 1.0, minimum purity 1.0.
- CUR's BRANCH clusters (allele-level SNPs only, est_div < delta) are never `copy` calls.

| bar | value | |
|---|---|---|
| H1 | ABSENT@0.02 12 / 12, BRANCH 12 / 12, ABSENT@0.01 8 against CUR's 8 | pass |
| H2 | 0 of 60 | pass |
| H3 | BRANCH +12; false-call families -29 | pass |
| H4a | TWIN_ALT, R0 9 -> S1 0 | pass |
| H4b | BRANCH, S1 0 -> S12 12 | pass |
| H4c | OUTSIDE, S12 20 -> S123 0 | pass |

#### Verdict (2026-10-15; corrected after the whole-branch review)

**Reading 1: claimed, in simulation.** H1-H4 pass as pre-registered. In ideal scenarios built from curated CAT isoforms, the
graph-informed O3 detector (S123) does three things:
- it finds absent copies that carry a novel exon whose sequence is absent from the family's loci and the genome: BRANCH
  12 / 12, against 0 / 12 for the current method;
- it stops false calls from uncatalogued paralogs (OUTSIDE 20 -> 0) and from copies a representative makes look identical
  (TWIN_ALT 9 -> 0);
- it keeps the current method's count of sequence-divergent absent copies: 12 / 12 at 2%, 8 / 12 at 1%.

**What each step does here.**
- Step 2 is needed for BRANCH and step 3 for OUTSIDE.
- Step 1 removes TWIN_ALT's false calls on its own (R0 9 -> S1 0, H4a), but step 3 removes them too. T is in the genome, so its
  reads' median `de` is only 0.001-0.002.
- The detector without step 1 (S23 = S3 plus the kept branch calls, an arm added after the review) has the same outcome as S123
  in every family of both the held-out and development sets. **So in these worlds step 1 is redundant with step 3.**
- H4a tests step 1 against R0, not whether S123 needs it. Step 1's own value (PSVs from every isoform, the MCL117 case) is for
  O2-style assignment and for cases step 3 cannot see; neither was tested here.

**Costs and limits.**
- **Step 3 rejected a correct 1% call on both sets.** Development SALL4 (median `de` 0.0095) was a detection lost.
  Held-out PPDPF was a pure 47-read call removed (median `de` 0.0059), though it missed the recall bar anyway. The reads' `de`
  and the cluster's est_div are different estimators held to the same delta; at 1% divergence they straddle it.
- **The 1% absent copies sit at delta**, and every arm misses the same 4 of 12:
  - three are pure clusters labelled allele-level (est_div 0.0069-0.0080);
  - the fourth (PPDPF) is a `copy` cluster holding 47 of the copy's 61 reads, under the 0.8 recall bar.
  - In S123 and S3, step 3 also removed that call.
- **The BRANCH scope is narrow.** The vg graph holds the whole loci, introns included (`vg rna -r` adds paths, it does not
  remove introns), and GraphAligner aligns to the graph, not to the paths. So a branch read means one with >= 50 bases that do
  not align to the catalogued loci. Structural changes made from sequence already present were not tested: an exon formed from
  intron sequence, a lost exon, a duplicated exon, a shifted splice site.
- **Recall and purity** (median 1.0) are over the absent copy's PLACED reads. In BRANCH_PPDPF, minimap2 left 30 of E's 61
  reads unmapped. The call holds 31 / 31 placed reads, 0.51 of all of E's reads.
- **Not all of CUR's BRANCH clusters are allele-level SNPs.** Two reach est_div 0.0073 (LINC01370, 59 E + 16 A reads) and
  0.0071 (ZBP1), likely misalignment around the novel exon. None is a `copy` call.
- **Read names encode the source copy** (famsim `copy|isoform|i`), and the clustering breaks coverage ties by name, in the
  prototype and in Rust alike. The reviewer's re-run with hashed names changed no outcome. Future worlds should anonymise
  names.
- **Scope:**
  - famsim reads from curated isoforms in small worlds; step 3's genome is the world's genome, not a whole genome.
  - Real data's other failure modes (other species, fragmented clusters) are not modelled.
  - "-29 false-call families" is a count over scenarios built to provoke those calls, not a rate.

**Changed after the results (no outcome changed on either set):**
- Step 3 now subtracts what a read's supplementary alignments cover from its clipped bases. A readthrough's clipped part
  aligned elsewhere is explained by the genome.
- Step 3 counts a call member without a primary alignment as wholly unexplained.
- The S23 arm was added.
- The detector's description is corrected: S1 aligns every isoform with minimap2 (the same isoforms the graph carries as
  paths), and step 2 asks whether bases align to the loci graph.


---

## PREREG_o3_graph_union2_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_union2_2026-10-15.md`)*

## Union-path graph O2 with padded path ends, and family-pooled O3, on a fourth handpicked family set (2026-10-15)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

Round 3. Round 2 (`docs/PREREG_o3_graph_union_2026-10-15.md`) claimed union O2 in part and graph O3 without the removed-copy
clause. Each failure pointed at one fix:
- **U4** (0.981 on paralog-structure reads): 5 AMY1 reads ran 16 bp past the end of union O2's own projection of their
  structure. Those bases counted as distance against the true copy, more than the near-identical copies' PSVs.
- **G3** (2 of 3 removed copies): PSG's removed copy was split over four hosts (27, 16, 10, 8 reads). Each host's cluster was
  pure, but none held 80% of the copy's reads.

### Methods under test

**O2, primary: padded union O2** (`o2graph.py run WORLD unionp 2`). Round 2's union O2 (every family isoform structure
threaded through every copy), with each path's two transcript ends extended by 50 bp of the copy's own genome. A read
running past an annotated end is compared with each copy's real sequence there. It is neither penalised nor excused for
the annotation.

A first version, which simply did not count read bases past a path end ("unione"), fixed AMY. But it also dropped real end
evidence: on dev-3 CROSS + REMOVED it fell from x1.07 to x1.03 and made 2 wrong calls. Padding replaced it.

**O3, primary: H** (`detect.py`). As round 2's G (cluster only the reads union O2 cannot explain, the mosaic rule, branch
calls, the genome check), but on padded union O2's unexplained reads POOLED over all of the family's hosts. The pool is
aligned to each host and clustered there, so an absent copy that copy_assign spreads over several hosts clusters once.

Reported beside them: the current O2 and CUR; round 2's union O2 and G; path O2, locus O2; S23, S123m; padded union O2's
gains split into PSV-separated (the runner-up copy can carry the read's structure) and structure-separated (`o2graph.py
gains`; round 2: 35% / 65%).

### Worlds

As round 2, seed 41: INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon), and CROSS with
`--projection liftover` (paralog structures lifted through the loci's genomic alignment, independent of union O2's own
spliced projection).

### Families (held out, picked before any run)

X-linked cancer-testis families SSX, MAGEA, SPANX; UGT2B; the opsins OPN1LW/MW; CGB; KIR; CYP4F. TRIM, LILR, DEFB10 and CSAG
were named as spares, in order, should one form no family; none was needed.

| family | copies / isoforms |
|---|---|
| SSX | 5 / 10 |
| MAGEA | 5 / 16 |
| SPANX | 4 / 4 |
| UGT2B | 3 / 8 |
| OPN1 | 3 / 9 (OPN1LW, OPN1MW, OPN1MW2) |
| CGB | 6 / 10 |
| KIR | 8 / 16 |
| CYP4F | 2 / 8 |

`/mnt/linuxdisk/tmp/o3_graph/families3.json` (its earlier sets are identical to `families2.json`).

### Known before running (development: dev, dev-2 = register 1298's held-out set, dev-3 = round 2's held-out set)

O2, padded union O2 (margin 2) against the current O2 and round 2's union O2; correct / assigned, wrong-family calls counted:

| set | worlds | current | **padded union O2** | union O2 (round 2) |
|---|---|---|---|---|
| dev | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced), CROSS (liftover) | 10,939 / 10,957 | **11,315 / 11,315 (x1.034)** | 12,391 / 12,391 (x1.133) |
| dev-2 | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced) | 11,530 / 11,657 | **12,701 / 12,701 (x1.102)** | 12,701 / 12,701 (x1.102) |
| dev-3 | INTACT, CROSS (liftover), OUTSIDE, REMOVED, BRANCH | 6,678 / 6,678 | **7,106 / 7,106 (x1.064)** | 7,170 / 7,175 (x1.074) |

- Paralog-structure reads, padded union O2: dev 883 / 883, dev-2 611 / 611, dev-3 261 / 261 (round 2's AMY1 errors gone).
- Padding costs most on dev (x1.133 -> x1.034). There copies' annotated isoforms end at different points, and the reads,
  trimmed inside their own isoform, ran past a paralog's shorter path end. That counted as a difference, but it is the
  annotation again, not sequence. Padding lets each copy's genome decide those bases, so the honest gain is family-set
  dependent (x1.03-x1.10) -- hence V2's 2% below (5% and 10% reported).
- Gain split, padded union O2 on dev-3: separated by PSVs 102, by structure 342 (the runner-up cannot carry the read's
  structure under the projection rules).

O3 (families with a detection or a false call; detectable = median `de` > delta):

| set | | CUR | S23 | G | **H** |
|---|---|---|---|---|---|
| dev | BRANCH detected | 0 / 8 | 8 / 8 | 8 / 8 | **8 / 8** |
| dev | REMOVED detected / detectable | 3 / 4 | 3 / 4 | 4 / 4 | **4 / 4** |
| dev | false-call families, INTACT + OUTSIDE + both CROSS | 4 / 32 | 0 | 0 | **0** |
| dev-2 | BRANCH | 0 / 8 | 8 / 8 | 8 / 8 | **8 / 8** |
| dev-2 | REMOVED | 0 / 2 | 1 / 2 | 2 / 2 | **2 / 2** |
| dev-2 | false-call families of 40 | 13 | 0 | 0 | **0** |
| dev-3 | BRANCH | 0 / 8 | 8 / 8 | 8 / 8 | **8 / 8** |
| dev-3 | REMOVED | 2 / 3 | 2 / 3 | 2 / 3 | **3 / 3** (PSG, pooled) |
| dev-3 | false-call families of 40 | 3 | 0 | 0 | **0** |

### Bars (held out: the 8 families above, seed 41, five worlds; wrong-family calls count as wrong)

| id | bar |
|---|---|
| V1 | padded union O2 precision >= 0.995 in every world (reads of catalogued copies) |
| V2 | padded union O2 correct >= 1.02 x the current O2's, pooled over the five worlds |
| V3 | padded union O2 correct >= the current O2's in every world |
| V4 | padded union O2 precision >= 0.99 on CROSS's paralog-structure reads |
| V5 | padded union O2 assigns no more reads of non-catalogued copies than the current O2, pooled |
| H1 | H detects >= 6 of the 8 BRANCH absent copies |
| H2 | H has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| H3 | H detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| H4 | H detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |

Detectable (H3): the removed copy's reads have a median primary `de` on the masked genome > delta. Detection and false
calls use `score.py`'s rule (register 1294).

### Readings (fixed now)

1. **V1-V5 pass:** claimed, "on simulated reads of handpicked real families, O2 on a family graph that threads every isoform
   structure through every copy, each path ending in its copy's own genome, assigns more reads correctly than the current
   O2 (>= 2%, the figure reported) at >= 99.5% precision, including when isoforms are annotated on another paralog".
2. **H1, H2, H4 pass and H3 passes:** claimed, "graph O3 that clusters the family's O2-unexplained reads together detects
   novel-exon absent copies and >= 75% of detectable removed copies, with at most one false call over intact,
   uncatalogued-paralog and paralog-annotation families". H3 not measurable: claimed without the removed-copy clause.
3. Whether the O2 step adds detections: H detects more than S23 (H4 holding) -> "the O2-informed pooling adds detections on
   held-out families"; equal -> "it adds nothing on this set" (reported, not a failure).
4. Each claim stands or falls separately. A failed bar is reported by name.

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome. Real data's other species,
fragments and allele background are not modelled.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (held out, 2026-10-15)

Run: code `c9a02337`, plus `ea2f4197` (below), the frozen copy_assign binary, seed 41. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout3/`.

**One change during the run, disclosed.** `vg msga` (v1.73) aborted (signal 6) on UGT2B's three loci, which the branch calls
need. The detector now retries a crashing `vg msga` with `-a` (build in FASTA order). That is the only option that built the
graph, and it is used only after a crash (`ea2f4197`, test). No other family, development or held out, crashed. The
detector was then re-run on INTACT, and the other four worlds ran with it.

**O2** (correct / assigned; there were no wrong-family calls):

| world | reads | current | **padded union O2** | union O2 (round 2) | path O2 | locus O2 |
|---|---|---|---|---|---|---|
| INTACT | 2,167 | 1,921 / 1,921 | 1,920 / 1,920 | 1,941 / 1,941 | 2,044 / 2,044 | 1,827 / 1,827 |
| CROSS (liftover) | 2,473 | 2,143 / 2,143 | 2,107 / 2,107 | 2,116 / 2,119 | 2,046 / 2,064 | 2,046 / 2,046 |
| OUTSIDE | 1,685 | 1,621 / 1,621 | 1,619 / 1,619 | 1,619 / 1,619 | 1,619 / 1,619 | 1,517 / 1,517 |
| REMOVED | 1,685 | 1,624 / 1,624 | 1,619 / 1,619 | 1,619 / 1,619 | 1,619 / 1,619 | 1,517 / 1,517 |
| BRANCH | 2,167 | 1,921 / 1,921 | 1,920 / 1,920 | 1,941 / 1,941 | 2,044 / 2,044 | 1,827 / 1,827 |
| **pooled** | 10,177 | 9,230 / 9,230 | **9,185 / 9,185 (x0.995)** | 9,236 / 9,239 | 9,372 / 9,390 | 8,734 / 8,734 |

- CROSS's paralog-structure reads (306): padded union O2 187 / 187, current 222 / 222, path O2 2 / 20, union 175 / 178.
- Reads of non-catalogued copies assigned, pooled: padded union O2 375, current 891.
- Counting those as the wrong calls they are, precision is 0.961 for padded union O2 and 0.912 for the current O2 (round 2's
  set: union O2 0.969, current 0.917).
- **No headroom.** 9,242 reads are separable by the genome at all (`seqsupport`), and the current O2 gets 9,230 of them right.
  - In CYP4F, KIR, OPN1 and UGT2B every method gets every read.
  - In MAGEA and SSX every method leaves about 40% unassigned (copies identical over the reads: MAGEA2 / MAGEA2B, SSX4 /
    SSX4B).
  - Padded union O2 loses 45 reads, called `unexplained`: 20 MAGEA and 15 SSX paralog-structure reads in CROSS, whose
    lifted structure no path matches within delta, plus SPANX and CGB reads.
- Gain split: 9 PSV-separated, 3 structure-separated.

**O3** (8 families per world):

| | CUR | S23 | S123m | G | **H** |
|---|---|---|---|---|---|
| BRANCH detected | 0 | 8 | 8 | 8 | **8** |
| INTACT + OUTSIDE + CROSS false-call families (of 24) | 2 | 0 | 0 | 0 | **0** |
| REMOVED detected (4 detectable) | 1 | 1 | 1 | 1 | **2** |
| REMOVED + BRANCH false-call families | 0 | 0 | 0 | 0 | **0** |

The removed copies' median `de` on the masked genome:
- detectable: UGT2B 0.074, CYP4F 0.047, KIR 0.024, CGB 0.0101;
- within delta: OPN1 0.0052, SPANX 0.0055, SSX 0.0016, MAGEA 0.0015.

H detects UGT2B, as every arm does, and CGB, which only H detects: copy_assign spreads CGB's 60 reads over two hosts (48 +
12), and the pool calls them once.

| bar | value | |
|---|---|---|
| V1 | precision 1.0 in every world | pass |
| V2 | 9,185 vs 9,230: x0.995 | **fail** |
| V3 | below the current O2 in every world (by 1-36 reads) | **fail** |
| V4 | 187 / 187 | pass |
| V5 | 375 vs 891 | pass |
| H1 | 8 / 8 | pass |
| H2 | 0 of 24 | pass |
| H3 | 2 of 4 detectable (0.5); CUR's one (UGT2B) detected | **fail** |
| H4 | H 10 vs S23 9, false-call families 0 = 0 | pass (H detects more: the O2-informed pooling adds detections) |

**What the misses show (post hoc).**
- KIR: copy_assign hosts all 60 reads on one copy, and union O2 calls all 60 `unexplained`. Yet the host cluster counts as
  explained: some catalogued copy carries every consensus allele at the cluster's sites. The `explains` rule checks only the
  sites where >= 3 reads differ from the host, so it misses the copy's other differences from the cluster (the reads are
  2.4% from the genome). Requiring the nearest single copy within delta before "explained" would call it.
- CYP4F (2 copies): no read of the removed copy is hosted on the family, so no arm sees it.

#### Verdict (2026-10-15)

*Superseded by "Verdict, revised after the review" below.*

**O2: not claimed (V2, V3 fail; V1, V4, V5 pass).** On this set the current O2 already reaches the genome-separable ceiling,
and padded union O2 assigns 0.5% fewer reads correctly. What holds on both fresh sets (rounds 2 and 3):
- graph O2 never assigns a read to a wrong copy or a wrong family (precision 1.0 / 0.9993);
- it stays right when isoforms are annotated on another paralog (187 / 187; path O2 2 / 20);
- it places 58-63% fewer reads of non-catalogued copies (375 vs 891, 223 vs 608), so its precision counting those is
  0.96-0.97 against 0.91-0.92.

Its gain in correct assignments is family-dependent: +6 to +10% where the current O2 leaves separable reads unassigned
(development sets), about 0 where it does not (this set).

**O3: claimed without the removed-copy clause (H1, H2, H4 pass; H3 fails).** Graph O3 replicates on its third fresh family
set:
- novel-exon absent copies 8 / 8 (current method 0);
- no false call over 24 families (current method 2);
- pooling the O2-unexplained reads across the family adds a removed-copy detection (CGB) that no per-host arm makes;
- it misses KIR (the `explains` rule, above) and CYP4F (no hosted reads).

### Review corrections (rounds 2-3 review, 2026-10-15; post hoc)

A fresh reviewer re-derived every number above from the raw outputs; all reproduce except the dev-2 false-call count (M1
below). The other corrections are about what the numbers mean.

**Verdict labels (C1).** Reading 2 allows "without the removed-copy clause" only when H3 is not measurable. H3 was measurable
(4 detectable) and failed. So: **O3 not claimed (H3 fails).** O2 was already not claimed.

**The KIR explanation was wrong (C2).** The `explains` rule did not explain KIR's cluster. Its status was `unexplained`: 5 of
41 alt sites copy-private, linkage 0.98, but est_div 0.0032, below delta, so no PSV call. The documented mosaic rule (this
prereg: "a cluster no single copy explains whose nearest copy differs from its consensus at > delta") would have called it,
since the nearest single copy is 0.012-0.032 away, depending on the host. But the code applied it only to `unresolved` / `unlinked` clusters. The
code now follows the documentation (`b0a3e59f`, test). Re-run with it, post hoc, which cannot change the verdict:
- here: H detects UGT2B, CGB and KIR, 3 of 4 detectable (0.75, H3's bar), with no false-call family in any world;
- on round 2's set: 3 of 3 (FCGR, LRRC37A, PSG), no false call;
- the frozen calls are in each held-out world's `graph/calls.frozen.tsv`.

**H flagged the same pooled reads once per host (I5).** For example, PSG's 61 reads appeared as 6 calls. The family-level
scorer counted no extra detection or false call, but "clusters once" was wrong. Calls are now one per member set
(`b0a3e59f`). Pooling also raises false calls before the genome check: Hn has false calls in CROSS (MAGEA, SSX) and OUTSIDE
(CGB, MAGEA, UGT2B; KIR too with the post-review detector) that only the genome check removes. H2 rests on the genome check.

**The O2 statements across rounds overclaimed (I1).**
- Precision: 1.0 on this set and 0.9993 for round 2's union O2. Round 2's 5 wrong AMY calls failed U4. Padded union's 1.0 on
  round 2's set is development data.
- Non-catalogued placements, 375 vs 891 (round 2: 223 vs 608), come mostly from the shared delta abstention. Locus O2, with
  no annotation and the same `decide()` rule, also places 48-55% fewer (467 vs 891; 273 vs 608).
- "187 / 187 on paralog-annotated isoforms" is right by abstaining: 119 of the 306 such reads are left unassigned.

**The liftover CROSS shares union's chains here too (I2).** 21 of the 26 lifted structures carry an intron chain union
threads. The 5 that differ carry 53 reads: union and padded union assign none of them (all `unexplained`), and the current O2
gets 35 right. Those reads are the Hn false calls in CROSS.

**Headroom wording (I4).**
- The share of padded union's net gain that the genome alone separates is 25% (dev), 85% (dev-2) and 23% (dev-3), so the
  sentence "where the current O2 leaves separable reads unassigned" fits only dev-2.
- Dev's padded gain is +3.4%.
- On this set, path O2 makes 372 correct calls the genome cannot separate (+142 net, mostly MAGEA), from the per-copy
  annotation, which union by design does not use.
- Union's structure-based separation found 3 distinct reads here, and 12 genome-separable reads remain unassigned.

**Smaller corrections.**
- **M1.** Dev-2's false-call count for CUR is 16 of 40, not 13 (13 was the 32-family count).
- **M2.** The 45 lost reads are net: 57 lost, all `unexplained` (MAGEA 22, SSX 15, CGB 10, SPANX 10), minus 12 gained (3
  distinct CGB reads). "Every method gets every read" is false for KIR (locus O2 1,900 of 2,336; path O2 2,285; union 2,326)
  and UGT2B (path O2 785 of 809; union 806 of 809), and CYP4F leaves 122 of 488 reads to every method (366 of 488). MAGEA and SSX's ~40% unassigned holds in INTACT, CROSS and BRANCH.
- **M3.** CYP4F (2 copies): in OUTSIDE and REMOVED copy_assign writes no row for the family (one catalogued copy left), so it
  was unreachable by design and still in H3's denominator. Round 4 requires >= 3 copies.
- **M4.** All five worlds share INTACT's catalogued-copy reads; pooled counts repeat reads up to five times.
- **M6.** CGB is at the boundary: median `de` 0.0101 against delta 0.00958, and it is detected only when the pool is aligned to
  host 1 (12 hosted reads), not host 0 (48).

#### Verdict, revised after the review (2026-10-15)

**O2: not claimed (V2, V3 fail).**
- No gain in correct assignments on this set.
- Precision 1.0.
- Fewer non-catalogued placements, mostly from the delta abstention it shares with locus O2. It keeps 5% more correct calls
  than locus O2.

**O3: not claimed (H3 fails, 2 of 4).**
- Branch calls: 8 / 8; no false call after the genome check.
- The documented mosaic rule, which the code skipped, would have given 3 of 4. That is post hoc; it is round 4's test.


---

## PREREG_o3_graph_union_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_union_2026-10-15.md`)*

## Union-path graph O2 and union-O2-informed O3 on a fresh set of handpicked families (2026-10-15)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

This follows `docs/PREREG_o3_graph_real_2026-10-15.md` (register 1298-1299). There, path O2 (each read against every curated
isoform of every copy) gained +14.3% only because each copy's catalog entry listed exactly the isoforms it expressed. When
copies also expressed isoforms annotated on a paralog (world CROSS), it named the annotated paralog: precision 0.42 on those
reads, 0.946 overall, and fewer correct calls than the current O2. On O3, the O2 step added nothing (S4 = S23).

### Methods under test

**O2** (`o2graph.py`; the reads O2 places on a family = copy_assign's rows):

| method | what it does |
|---|---|
| current O2 | copy_assign `--no-as-tied-only`: rows with status `assigned` |
| **union O2** | the family graph threads EVERY family isoform structure through EVERY copy. A copy's paths: its own curated isoforms, plus each paralog isoform spliced-aligned to its locus (kept with >= 95% coverage, the same exon count, canonical junctions, an intron chain new to it). Then as path O2: edits + unaligned bases to the copy's best path; assigned when the next copy is >= MARGIN farther, abstains when tied or thin, `unexplained` beyond delta. **MARGIN = 2.** Only the sequence separates copies, not which copy the annotation lists a structure on. |
| path O2 | each copy's own curated isoforms only (register 1298), reported |
| locus O2 | the annotation-free baseline: the same decision on spliced alignments to each copy's genomic locus (any canonical junction), reported |

**O3** (`detect.py`):

| arm | what it does |
|---|---|
| CUR | `copy_assign --psv-clusters` |
| S23 | branch calls (reads on no path of the family's vg graph) + PSV calls, kept only when the genome does not explain their reads |
| S123m | S123 with the MOSAIC rule: a cluster no single copy explains whose NEAREST copy differs from its consensus at > delta of its bases (GOLGA8's removed copy, register 1298, carries alleles of several copies) |
| **G** | union-O2-informed O3: the clustering run only on the reads union O2 calls `unexplained`, with the mosaic rule, plus the branch calls, kept by the genome check |

### Worlds (`realworld.py`)

As register 1298 (handpicked CAT families on the real CHM13 genome, reads from curated isoforms, anonymised names), seed 31:
INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon, 0.2% SNPs), and **CROSS with
`--projection liftover`**. In this CROSS, every copy also expresses up to 2 paralog-annotated structures, lifted onto it through
the two loci's genomic alignment (minimap2 asm20). That projection is independent of the spliced transcript alignment union O2
itself uses. On development, 45 of its 66 structures matched the spliced projection exactly, 12 differed in exon ends, and 9
were different choices.

### Families (held out, picked before any run)

Soto's human-specific duplications FAM72, ARHGAP11 and TCAF2; copy-number mainstays AMY, SULT1A, FCGR, APOBEC3 and PSG.
ARHGAP11 and APOBEC3 form no family at the build's 90% identity linkage (built, no world). GTF2H2 (the SMA region) and LRRC37A
(17q21.31) were added in their place, still before any run.

| family | copies / isoforms |
|---|---|
| FAM72 | 4 / 11 |
| TCAF2 | 2 / 6 (TCAF2, TCAF2P1) |
| AMY | 9 / 30 |
| SULT1A | 4 / 16 |
| FCGR | 4 / 16 |
| PSG | 7 / 28 |
| GTF2H2 | 2 / 8 |
| LRRC37A | 3 / 8 |

`/mnt/linuxdisk/tmp/o3_graph/families2.json` (its dev and first held-out sets are identical to `families.json`).

### Known before running (development)

Development = the 8 original development families, plus the 8 families of register 1298's held-out set. Those are seen now
and count as development ("dev-2").

O2, union O2 (margin 2) against the current O2; correct / assigned, wrong-family calls counted:

| set | worlds | current | union O2 | path O2 | locus O2 |
|---|---|---|---|---|---|
| dev | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced) | 8,703 / 8,712 | **9,944 / 9,944 (x1.143)** | 9,967 / 10,203 | 8,612 / 8,612 |
| dev | CROSS (liftover) | 2,236 / 2,245 | **2,447 / 2,447 (x1.094)** | 2,222 / 2,443 | 2,247 / 2,247 |
| dev-2 | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced) | 11,530 / 11,657 | **12,701 / 12,701 (x1.102)** | 12,707 / 12,865 | 11,285 / 11,436 |

- Paralog-structure reads: union O2 is right on every one it assigns (dev 490 / 490, liftover 395 / 395, dev-2 611 / 611).
  Path O2 gets 67 / 303, 66 / 287 and 115 / 273.
- Locus O2 gains nothing (x0.99, x0.98). The gain is not genome alignment as such but the family's structures with their PSVs.
- Over dev-2's four original worlds alone union O2 is +9.7%, below register 1298's 10% bar; hence U2's 5% below.

O3, families with a detection or a false call:

| set | scenario | CUR | S23 | S123m | G |
|---|---|---|---|---|---|
| dev | BRANCH detected | 0 | 8 | 8 | 8 |
| dev | REMOVED detected (4 detectable) | 3 | 3 | 3 | **4** (GAGE: the O2 restriction) |
| dev | false-call families, INTACT + OUTSIDE + CROSS | 4 | 0 | 0 | 0 |
| dev-2 | BRANCH detected | 0 | 8 | 8 | 8 |
| dev-2 | REMOVED detected (2 detectable) | 0 | 1 | 2 | **2** (GOLGA8: the mosaic rule) |
| dev-2 | false-call families of 32 | 13 | 0 | 0 | 0 |

### Bars (held out: the 8 families above, seed 31, five worlds; wrong-family calls count as wrong)

| id | bar |
|---|---|
| U1 | union O2 precision >= 0.995 in every world (reads of catalogued copies) |
| U2 | union O2 correct >= 1.05 x the current O2's, pooled over the five worlds |
| U3 | union O2 correct >= the current O2's in every world |
| U4 | union O2 precision >= 0.99 on CROSS's paralog-structure reads (register 1299's X1, which path O2 failed) |
| U5 | union O2 assigns no more reads of non-catalogued copies than the current O2, pooled |
| G1 | G detects >= 6 of the 8 BRANCH absent copies |
| G2 | G has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| G3 | G detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| G4 | G detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |

Detectable (G3) is as before: the removed copy's reads have a median primary `de` on the masked genome > delta. Detection and
false calls use `score.py`'s rule (register 1294).

### Readings (fixed now)

1. **U1-U5 pass:** claimed, "on simulated reads of handpicked real families, O2 on a family graph that threads every isoform
   structure through every copy assigns >= 5% more reads correctly than the current O2, at >= 99.5% precision, including
   when isoforms are annotated on another paralog".
2. **G1, G2, G4 pass and G3 passes:** claimed, "graph O3 detects novel-exon absent copies and >= 75% of detectable removed
   copies, with at most one false call over intact, uncatalogued-paralog and paralog-annotation families". G3 not measurable:
   claimed without the removed-copy clause.
3. Whether the O2 restriction and the mosaic rule ADD detections: G detects more than S23 (G4 holding) -> "they add
   detections on held-out families"; equal -> "they add nothing on this set" (reported, not a failure).
4. Each claim stands or falls separately. A failed bar is reported by name.

Reported beside the bars: path O2 and locus O2; G against S123m, S123 and S23 per family; union O2 against register 1298's 10%;
union O2's gains split into reads whose runner-up copy also carries the read's structure (separated by PSVs) and reads whose
runner-up copy cannot carry it (separated by structure).

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome. Real data's other species,
fragments and allele background are not modelled.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (held out, 2026-10-15)

Run: code and prereg `4850047a`, the frozen copy_assign binary, seed 31. Outputs in `/mnt/linuxdisk/tmp/o3_graph/real/heldout2/`
(one directory per world, `o3_score.json`, `o2sum_*.json`).

**O2** (correct / assigned; wrong-family calls counted; there were none):

| world | reads | current | **union O2** | path O2 | locus O2 |
|---|---|---|---|---|---|
| INTACT | 2,125 | 1,452 / 1,452 | 1,561 / 1,561 (x1.075) | 1,650 / 1,650 | 1,473 / 1,473 |
| CROSS (liftover) | 2,635 | 1,706 / 1,706 | 1,816 / 1,821 (x1.064) | 1,662 / 1,900 | 1,736 / 1,736 |
| OUTSIDE | 1,640 | 1,034 / 1,034 | 1,116 / 1,116 (x1.079) | 1,150 / 1,150 | 1,032 / 1,032 |
| REMOVED | 1,640 | 1,034 / 1,034 | 1,116 / 1,116 (x1.079) | 1,150 / 1,150 | 1,032 / 1,032 |
| BRANCH | 2,125 | 1,452 / 1,452 | 1,561 / 1,561 (x1.075) | 1,650 / 1,650 | 1,473 / 1,473 |
| **pooled** | 10,165 | 6,678 / 6,678 | **7,170 / 7,175 (x1.074, precision 0.9993)** | 7,262 / 7,500 (0.968) | 6,746 / 6,746 |

- CROSS's paralog-structure reads (510): union O2 255 / 260 correct (0.981), current 254 / 254, path O2 12 / 250.
- Reads of non-catalogued copies assigned, pooled: union O2 223, current 608, path O2 250.
- Per family (correct, current -> union): FCGR 511 -> 851, SULT1A 772 -> 867, AMY 646 -> 714, GTF2H2 357 -> 366,
  LRRC37A 810 -> 810, TCAF2 384 -> 384, FAM72 1,077 -> 1,070, PSG 2,121 -> 2,108.
- Genome-separated reads (6,775): current 6,666, union O2 6,760, path O2 6,508 correct.

**O3** (8 families per world):

| | CUR | S23 | S123m | **G** |
|---|---|---|---|---|
| BRANCH detected | 0 | 8 | 8 | **8** |
| INTACT + OUTSIDE + CROSS false-call families (of 24) | 3 | 0 | 0 | **0** |
| REMOVED detected (3 detectable) | 2 | 2 | 2 | **2** |
| REMOVED + BRANCH false-call families | 0 | 0 | 0 | 0 |

The removed copies' median `de` on the masked genome: FCGR 0.053, PSG 0.042 and LRRC37A 0.021 are detectable. TCAF2 0.0070,
GTF2H2 0.0062, FAM72 0.0035, AMY 0.0017 and SULT1A 0.0015 are within delta. Every arm detects FCGR and LRRC37A and misses PSG.

| bar | value | |
|---|---|---|
| U1 | lowest world precision 0.9973 (CROSS) | pass |
| U2 | 7,170 vs 6,678: x1.074 | pass |
| U3 | union >= current in all five worlds | pass |
| U4 | 255 / 260 = 0.981 | **fail** |
| U5 | 223 vs 608 | pass |
| G1 | 8 / 8 | pass |
| G2 | 0 of 24 | pass |
| G3 | 2 of 3 detectable (0.67) | **fail** |
| G4 | G 10 = S23 10, false-call families 0 = 0 | pass (the O2 restriction and the mosaic rule add nothing on this set) |

**Gain split (pre-registered as reported).** `o2graph.py gains`: union O2's correct calls where the current O2 abstains,
by whether the runner-up copy can carry the read's structure (its projection onto the runner-up's locus succeeds):

| world | separated by PSVs (runner-up carries the structure) | separated by structure (runner-up cannot) |
|---|---|---|
| INTACT | 42 | 68 |
| CROSS | 60 | 68 |
| OUTSIDE | 14 | 69 |
| REMOVED | 14 | 69 |
| BRANCH | 42 | 68 |
| total | 172 (35%) | 342 (65%) |

"Cannot carry" uses union O2's projection rules (>= 95% coverage, the same exon count, canonical junctions): the runner-up's
genome lacks an exon or a canonical splice site for that isoform. Locus O2, which allows any canonical junction on each
copy's genome, does not separate these reads (x1.010). So two thirds of the gain rests on the projection rules, i.e. on
the assumption that a copy without the structure's exons or splice sites does not make the isoform.

**What the failures show (post hoc).**
- U4: all 5 wrong calls are AMY reads (AMY.2 and AMY.5, two near-identical AMY1 copies) of one structure lifted from AMY.1.
  The lifted structure's last exon ends 16 bp after union O2's own spliced projection of the same isoform. The read's 16
  bases past the path end count as distance against the true copy, more than the copies' few PSVs. Read ends past a path
  end are annotation, not sequence: a decision over the span both copies align (end clips not counted) would remove it.
- G3: PSG's removed copy is detected in pieces. copy_assign hosts its 61 reads on four copies (27, 16, 10, 8). Each host's
  cluster is pure (est_div 0.021-0.034, reads' `de` 0.04 from the genome), but none holds 80% of the copy's reads. Union O2
  calls all 61 `unexplained`, yet G clusters them per host. Pooling the union-O2-unexplained reads across the whole family
  would call it once; register 1298's TBC1D3 calls (17 + 42 reads over two hosts) are the same split.

#### Verdict (2026-10-15)

*Superseded by "Verdict, revised after the review" below.*

**O2: claimed in part (U1, U2, U3, U5 pass; U4 fails).** On a fresh set of handpicked families:
- union O2 assigns 7.4% more reads correctly than the current O2, at precision 0.9993;
- it does better than the current O2 in every world, including when isoforms are annotated on another paralog (CROSS,
  structures lifted independently of its own projection);
- path O2 drops to 0.875 there;
- the paralog-structure clause fails narrowly (0.981 vs 0.99) on an end effect between two near-identical AMY1 copies;
- the gain is not plain genome alignment (locus O2 x1.010).

**O3: claimed without the removed-copy clause (G1, G2, G4 pass; G3 fails).** Graph O3 replicates on a third family set:
- novel-exon absent copies 8 / 8 (current method 0);
- no false call over 24 intact, uncatalogued-paralog and paralog-annotation families (current method 3);
- removed copies 2 of 3 detectable, the same as the current method;
- the O2 restriction and the mosaic rule add nothing here (G = S23).

**Next:** both failures point at one fix each, to be pre-registered on a new family set:
- decide over the span both copies align;
- pool the union-O2-unexplained reads across the family before clustering.

### Review corrections (rounds 2-3 review, 2026-10-15; post hoc)

A fresh reviewer re-derived every number above from the raw outputs, and all of them reproduce. The corrections are about
what the numbers mean.

**Verdict labels (C1).** The readings allow no partial claim. Reading 1 needs all of U1-U5. Reading 2 allows "without the
removed-copy clause" only when G3 is not measurable. U4 failed, and G3 was measurable (3 detectable) and failed. So:
**O2 not claimed (U4 fails); O3 not claimed (G3 fails).** The bars that passed are results, not claims.

**The O2 gain rests on one family, FCGR (I3).**
- 340 of the 342 structure-separated gains are FCGR: 68 distinct reads, counted in all five worlds.
- The runner-up copy, FCGR2C, has an AT-AG last intron. Union's canonical-junction rule rejects that projection, although it
  aligns 100% at NM 1. The read is 1 edit apart on the two genomic loci but 4 apart on union's paths.
- The simulator applies the same rule to CROSS's structures, so that separation is correct by construction.
- Without FCGR, union O2 is x1.025, short of U2's x1.05; padded union is x1.014.

**The gain split's PSV class is mislabelled (I3).** 70 of its 172 reads are AMY reads the genome does not separate. Annotated
end positions separate them, and they disappear with padding (round 3). The PSV-separated gain is 102 (round 3's dev-3
split).

**The liftover CROSS is mostly not independent of union's projection (I2).**
- 51 of its 52 lifted structures carry an intron chain union threads: 33 identical, 18 differing only in end position, by
  <= 17 bp.
- The one that differs (PSG, 9 reads) is called `unexplained` on all 9.
- So this CROSS tested tolerance to end positions, not to an independently projected structure.

**Overlapping worlds (M4).** All five worlds share INTACT's reads of catalogued copies. Pooled counts repeat the same reads up
to five times, so "every world" amounts to about three independent comparisons.

**Provenance (M5).** The gain split was computed with `o2graph.py gains`, added in `c9a02337` after this run's code
(`4850047a`). It re-aligns the run's saved read and path files.

**Boundary (M6).** LRRC37A's removed-copy detection holds exactly 48 of 60 reads, which is the recall bar (0.80).

**Held-out files.** `graph/calls.frozen.tsv` in each `heldout2` world is a copy of the 22:02-22:08 re-run (round-3 development), not of the 21:51 run; it matches all 440 rows these results were scored on and also carries H rows.
`graph/calls.tsv` now holds the post-review detector's calls.

#### Verdict, revised after the review (2026-10-15)

**O2: not claimed (U4 fails).**
- Union O2 assigned 7.4% more reads correctly at precision 0.9993 and passed U1, U2, U3 and U5.
- FCGR's single non-canonical intron carries that gain; without FCGR it is x1.025.
- Its paralog-annotation robustness was tested only against structures it threads itself.

**O3: not claimed (G3 fails, 2 of 3 detectable).** Branch calls and the genome check gave 8 / 8 novel-exon copies and no false
call over 24 families; those stand as results.


---

## PREREG_o3_psv_clusters_v2_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_psv_clusters_v2_2026-10-15.md`)*

## PREREG — O3 PSV clustering v2 (linkage test) in Rust (`copy_assign --psv-clusters`), held-out on base and hard worlds (2026-10-15)

**Written after development and before any held-out world below is built.** Follows `docs/PREREG_o3_unexplained_2026-10-15.md`
(v1, PASSED). Design: `docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md`.

### What changed since v1, and why (development only: seeds 32, 33, 51, 52)

- **Harder worlds** (`world.py --design hard`): the heterozygous allele background (H = A + 0.1 / 0.2 / 0.4%), a copy plus a background
  allele on one host, 8-copy families (intact and with one copy removed), gene conversion (an absent copy with two exons converted back
  from its host; present copies with converted exons), a mosaic absent copy, 5x sequencing error (a pessimistic stress case; Kinnex
  HiFi is ~0.1%), expression imbalance (host 200 reads, absent copy 6), exon-skipping isoforms, an absent copy that lost an exon,
  A-to-I RNA editing (8 sites, each in 40% of a present copy's reads), an absent copy at 4 reads, and an absent copy whose reads split
  between two hosts 0.4% apart. A 10x-error scenario was built and dropped (1% error is not a Kinnex regime).
- **v1 broke on editing** (false unexplained clusters of a present copy's reads: 3-9 reads, A->G at every site) **and on 1% error**
  (one 5-read cluster). **v2 adds the linkage test**: a missing sequence is a haplotype, so its alts travel together; edits and errors
  hit sites independently. Over every pair of a cluster's unexplained sites, P(alt at both | alt at either) among the host's reads
  covering both; median >= 0.8 or the cluster is `unlinked` (editing signature reported as the A->G / T->C fraction). Dev: every true
  missing sequence 0.95-1.00, editing 0.27-0.40, the error cluster 0.20.
- **Dev result of v2** (arm P, 4 worlds x E0/E2): absent copies >= 0.6% 207/208; control families 0/80; spurious unexplained clusters
  0; multi-group separation 16/16; alleles labelled `allele_or_copy` 70/70; copies >= 2% labelled `copy` 120/124 (the 4 are CONV_ABSENT,
  whose effective divergence is lowered by the converted exons). v1's base worlds unchanged (35/35, 43/43).
- **Rust port** (`rustle::missing_copy::psv_clusters`, `copy_assign --psv-clusters`; a post-pass over the run's own
  `<out>.assignments.tsv`, read sequences from one BAM sweep in their original orientation). **P0 parity: 1,241 / 1,241 clusters
  identical** to the prototype (membership by source, status, site counts, explaining copies, est_div / linkage / editing within
  rounding) on all 8 dev runs; the assignments it reads are byte-identical to the prototype's input; every other `copy_assign`
  product is byte-identical with and without the flag. 4 unit tests (cs parsing, Counter tie rule, a missing haplotype, unlinked edits
  and an explaining copy).

Frozen (md5): `psvclust.py` 17888625c3d31245d03e2221e5086dc0 (the reference), `score.py` 51b602dab31809c3743fb4d7aaacbcb8,
`world.py` b9b9cc581cdf845eb65001ed748a94d9, `parity.py` 858a8c4282b22fb2278a71831ba31b86; `copy_assign` sha1 f36d4e7a.

### Held-out run (the RUST product is scored: `score.py rows ... --rust`)

Base design seeds **61, 62, 63** and hard design seeds **71, 72, 73, 74, 75** (`--strands`), reads E0 and E2, arm P
(`copy_assign --families --flag-missing-copies --no-as-tied-only --psv-clusters`). Parity with the prototype is re-checked on every
held-out run (`parity.py`) and reported.

### Bars (pooled over all held-out runs)

| id | measure | bar | dev v2 |
|---|---|---|---|
| **H1** | absent copies with d >= 0.6% recovered (pure unexplained cluster) | >= 0.90 | 207/208 |
| **H2** | control families with an unexplained cluster (base: NULL_TWIN, NULL_MULTI, NULL_CLOSE, ADDBACK, DECOY; hard: MANY_NULL, NULL_CONV, HIGHERR_NULL, EDIT_NULL) — 152 units | <= 3 | 0/80 |
| **H2b** | spurious unexplained clusters (majority of a PRESENT copy's reads), all families | <= 3 | 0 |
| **H3** | families with two or three missing sequences on one host (COMBO_EH, DOUBLE, HET_COPY): all recovered separately | >= 0.80 | 16/16 |
| **H4** | recovered alleles labelled `allele_or_copy`; recovered copies with d >= 2% labelled `copy`, excluding CONV_ABSENT | >= 0.90 each | 70/70; 116/116 |
| **P0'** | held-out Rust-vs-prototype parity | every cluster identical | 1,241/1,241 |

### Reported beside (no bar)

Every scenario; the editing families (EDIT_NULL / EDIT_COPY: their unlinked clusters, linkage and A->G fraction); the 5x-error
scenarios; copies at 0.3% and 1%; CONV_ABSENT's labels; the existing detector (`--flag-missing-copies`) as the baseline.

### Reading, fixed now

- **H1-H4 and P0' pass** → `copy_assign --psv-clusters` (with `--no-as-tied-only`) is O3's method at simulation level, robust to the
  allele background, editing, conversion, many-copy families, imbalance and isoforms. Next level: real families on the main machine
  (remove copies from the GENOME, re-map, add back) with a single-copy-gene control for the real allele background, pre-registered.
- **P0' fails** → the Rust port has a defect: fixed in the port, never in the rule; the bars are read only after P0' passes.
- **H1, H2 or H2b fails** → v2 is not adopted; the failing scenarios are named; no constant is re-tuned on held-out worlds.
- **H3 or H4 fails alone** → detection is claimed, separation / labelling is not.

## OUTCOME (2026-10-15) — all bars PASS; the Rust product scored

`pixi run bash bench/o3_unexplained/heldout_v2.sh ~/rustle_data/o3_unexplained/heldout_v2` (base 61-63, hard 71-75, E0 / E2, arm P;
`copy_assign` sha1 f36d4e7a; frozen md5s as above).

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| P0' | Rust vs prototype, every held-out run | **16 / 16 runs, 0 cluster differences**; assignments byte-identical 16/16 | identical | ✅ PASS |
| H1 | absent copies >= 0.6% recovered | **425 / 426 = 0.998** (existing detector 414/426) | >= 0.90 | ✅ PASS |
| H2 | control families with an unexplained cluster | **0 / 152** | <= 3 | ✅ PASS |
| H2b | spurious unexplained clusters of a present copy's reads | **0** | <= 3 | ✅ PASS |
| H3 | two or three missing sequences on one host, all separated | **41 / 44 = 0.932** | >= 0.80 | ✅ PASS |
| H4 | alleles labelled `allele_or_copy`; copies >= 2% labelled `copy` (CONV_ABSENT excluded) | **135 / 135; 262 / 262** | >= 0.90 each | ✅ PASS |

**Reported beside.**
- RNA editing: EDIT_NULL 0/20 false flags, EDIT_COPY's absent copy 20/20; the editing produced 68 `unlinked` clusters (11 in
  EDIT_NULL, 57 in EDIT_COPY), none unexplained — the linkage test is what keeps them out.
- The one missed absent copy: HIGHERR@0.01 (5x error, a stress case), hard73 E2 f26 (its cluster was `unlinked`). HIGHERR@0.006
  20/20, labelled `allele_or_copy`.
- The 3 H3 misses are HET_COPY's background allele at 0.2% (17/20 recovered); its absent copy 20/20.
- The allele background: HET_BG at 0.1% recovered 14/30 (one or two variant sites: below what a cluster can call), 0.2% 28/30, 0.4%
  30/30, all labelled `allele_or_copy`; never a false `copy`.
- CONV_ABSENT labelled `copy` 20/20 here (4/8 `allele_or_copy` on dev: the label follows the conversion-reduced divergence).
- Labels at delta: copies at 1% split (TWIN 8/4, UNIQ 8/4, MULTI 8/4, MANY_REMOVED 12/8, DOUBLE 6/6, SPLIT 8/12 copy / allele_or_copy).
- Copies at 0.3% recovered 24/24 (existing detector 3/24).

**Reading (as fixed).** `copy_assign --families --no-as-tied-only --psv-clusters` is O3's method at simulation level: every missing
copy from 0.3% upward is recovered (425/426 at >= 0.6%), no intact family is flagged across the allele background, editing,
conversion, 8-copy families, imbalance, isoforms and 5x error (0/152), two or three missing sequences on one host are separated, and
the label states which ones RNA cannot tell from an allele. Next level: real families on the main machine (remove copies from the
GENOME, re-map the real reads, add back) with a single-copy-gene control for the real allele background, pre-registered before it runs.


---

## PREREG_o3_realdata_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_realdata_2026-10-15.md`)*

## O3 on real reads: copies removed from the genome, then added back (2026-10-15)

Pre-registration, written before any input of this test was built or run. Code: `bench/o3_realdata/` at `dd37275e`
(15 tests). Branch `o3-realdata`. Method under test: `copy_assign --psv-clusters` v2 (merged into `main` as `a493a66a`,
`docs/PREREG_o3_psv_clusters_v2_2026-10-15.md`), unchanged.

### Question

The advisor's O3 (2026-10-15): cluster the reads O2 places on a family by their PSVs; a cluster that no copy in the genome
explains better is a copy absent from the genome. In simulation it recovered 425 / 426 absent copies and flagged 0 / 152 intact
families. On real gorilla reads:

1. When a copy is removed from the GENOME and its reads are re-mapped, do they form an `unexplained` cluster on a surviving
   paralog?
2. Does the cluster go away when the copy is back in the genome?
3. Do single-copy genes, with the real heterozygous allele background, stay free of false `copy` calls?

This is step 6 of `docs/HANDOFF_laptop_to_main_2026-10-15.md`.

### Known before running (recorded now)

- **Simulation** (laptop, register 1292-1293): recovery 425 / 426 at ≥ 0.6% divergence, 0 / 152 control families flagged, RNA
  editing kept out. With exact AS ties only, the method sees twins only (80 / 310), so `--no-as-tied-only` is required.
- **The development substrate** (A13, `/mnt/linuxdisk/tmp/rna_allele/`) was built and used in 2026-10 for earlier O3 methods
  (IsoCon chain, `o3_candidates`; Amendments 7-15 of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md`):
  - 53 gorilla families of an older catalog (GWFAM), each with one copy hard-masked to N in the genome;
  - 59,013 real reads (up to 500 per copy) re-mapped to the masked genome (`linktest/R.bam`) and to the unmasked one
    (`control/R0.bam`);
  - the copy of origin of each read (`linktest/labels.tsv`).

  PSV clustering has never been run on it. Seen before: those methods' results (A13 detection 25 / 53 families; A14's control
  35 / 53 families with false flags), the families and their read counts.
- **Not seen**: any `--psv-clusters` output on real reads; which catalog families the held-out rule picks.
- **Allele cut-off**: delta = 0.00958, the 99th percentile of exonic divergence between KB3781's two haplotypes over its
  single-copy genes. So about 1% of single-copy genes carry alleles more divergent than delta. If expressed and clustered, such
  an allele is labelled `copy` by construction: that is the method's known limit, and the control measures its rate.

### Method (fixed)

**Removal and re-mapping.**
- The masked span is a copy's locus hull, hard-masked to N (`bedtools maskfasta`) in `GGO.fasta`.
- The reads are re-mapped with the baseline's exact command:
  `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (A13's `R.bam` header), sorted and indexed.
- The masked genome gets its own index (`minimap2 -x splice:hq -d`). The added-back arm maps the same reads to the unmasked
  genome with `winloci_data/GGO.splice.mmi`.

**The O3 run** (every arm):
`copy_assign --bam BAM --fasta GENOME --regions regions.txt --families COPIES.tsv --copies-fa COPIES.fa --no-as-tied-only
--psv-clusters --region-threads 1 --skip-poa-diagnostic --out PREFIX`.
The catalog is the surviving copies in the removed arm, and every copy in the added-back arm.

**Held-out families** (`prepare.py heldout`, seed 1):
- from the shipped KB3781 catalog (`runs/gorilla_KB3781/gorilla_KB3781.fam.copies.tsv`), families with ≥ 3 copies;
- every copy has ≥ 20 primary reads with an aligned base on its exons (a read on the exons of two copies counts for neither);
- every locus hull is ≤ 200 kb;
- no copy overlaps a copy of A13's 201;
- the copy last by (chrom, start) is masked, and a family whose masked hull overlaps any other catalog copy, its own siblings
  included, is dropped;
- if more than 60 families remain, a seeded sample of 60 is taken.

Reads: up to 500 per copy (seeded), labelled D (the masked copy) or S (a survivor), written in their original orientation
for `-uf`.

**Single-copy control** (`prepare.py control`, seed 1):
- candidates are KB3781's `S_single` genes with ≥ 500 exonic bp lifted ≥ 90% to the other haplotype; heterozygosity =
  (mismatches + indels) / lifted bp;
- the R3 set is a random sample of 300 expressed genes (≥ 20 primary reads on the exons), taken in seeded random order;
- also reported: up to 100 expressed genes with heterozygosity > delta (`high_het`), where `copy` calls are expected by
  construction;
- each gene is a one-copy family (its exons spliced from the genome), run on the baseline BAM and the unmasked genome.

**Scoring** (`score.py`; the simulation's rules, `bench/o3_unexplained/score.py`):
- A read is *placed* on a family when it has an assignment row of that family.
- **Recovered**: the `unexplained` cluster holding the most of the masked copy's placed reads holds ≥ 0.8 of them, with purity
  ≥ 0.9. Its label is `copy` when est_div > delta, else `allele_or_copy`.
- **Cleared** (added-back arm): no `unexplained` cluster of the family holds ≥ 0.5 of the masked copy's placed reads.
- **False copy** (control): an `unexplained` cluster with est_div > delta in a one-copy family.
- A masked copy is eligible when ≥ 6 of its reads are placed on its family in the removed arm.

### Bars (held-out families and the control; species: gorilla only)

| id | bar |
|---|---|
| R1 | recovered in ≥ 0.70 of the eligible masked copies |
| R2 | cleared in ≥ 0.90 of the recovered families |
| R3 | false `copy` in ≤ 0.02 of the 300 random single-copy genes |

**Development (A13's 53 families):** the same runs and scoring, reported without a bar. Its purpose is to exercise the
pipeline on a substrate that was built for this kind of test. If a defect shows up there, the fix is a recorded amendment
made before the held-out runs.

### Readings (fixed now)

1. **R1, R2 and R3 pass**: claimed, "O3 (PSV clusters) detects copies removed from the genome on real gorilla reads, the
   detection goes away when the copy is back, and single-copy genes stay below 2% false `copy` calls". The claim is quoted with
   the eligible and recovered counts.
2. **R1 passes, R2 or R3 fails**: detection reported; specificity not established, and the failing side is named.
3. **R1 fails**: not claimed; what happens to the masked copy's reads is reported (absorbed by a paralog, placed on another
   family, not placed), with the clusters they form.

Reported, not barred:
- the labels of the recovered clusters and the divergence of the masked copy from its nearest survivor;
- the masked copy's reads by fate;
- unexplained clusters per family in both arms;
- the control by heterozygosity bin, and the `high_het` stratum;
- the development substrate's numbers.

### Amendment 1 (2026-10-15, before any run; user decision)

Building the held-out inputs (`prepare.py heldout`, the rule above) picked only 14 families. Of the catalog's 90 families with
>= 3 copies, 62 fail the rule that every copy has >= 20 reads, 7 overlap A13 and 3 have a masked hull over another copy.
That leaves R1 coarse: 10 of 14 would pass. Seen at this point: the counts of families each rule removes, and nothing else.

**Change (the user chose it from three options):** the held-out families are those with **>= 2 copies** under the same
rules (`prepare.py heldout --min-copies 2`). 74 families pass, and a seeded sample of 60 is taken (seed 1). In a 2-copy
family, the masked copy's reads have one surviving paralog.

R1-R3 and the readings are unchanged. R1 and R2 are also reported separately for 2-copy and >= 3-copy families
(`score.py verdict --picked`). The 14-family build stays on disk (`/mnt/linuxdisk/tmp/o3_realdata/heldout/`) and is not run.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-15)

Run: code `39e333c1` (Amendment 1), `copy_assign` from `main` `a493a66a`. Outputs are in `/mnt/linuxdisk/tmp/o3_realdata/`
(`heldout2/`, `dev/`, `control/`).

**Cost.** The masked genome's index took 195 s at 21.2 GB. The 24,920 held-out reads were mapped in 4 parts per arm (about
220 s per part; one call of all of them hit the 10-minute limit). Each `copy_assign --psv-clusters` run took 166-478 s at
≤ 2.4 GB.

**Held out (Amendment 1's 60 families: 48 with 2 copies, 12 with >= 3).**
- Reads: 9,893 masked-copy (D) and 15,027 survivor (S). Masked arm: 24,904 primary and 16 unmapped; added-back arm: 24,920
  primary.

| | held out | development (A13, 53 families) |
|---|---|---|
| eligible masked copies (>= 6 reads placed on the family) | 12 (all >= 3-copy) | 25 |
| **R1: recovered** | **4 / 12 = 0.333, FAIL** (bar 0.70) | 15 / 25 = 0.600 |
| **R2: cleared after the copy is added back** | **3 / 4 = 0.750, FAIL** (bar 0.90) | 9 / 15 = 0.600 |
| labels of the recovered clusters | `copy` 3, `allele_or_copy` 1 | `copy` 10, `allele_or_copy` 5 |
| unexplained clusters, removed / added-back arm | 58 / 138 | 177 / 210 |

**One-copy families get no placed reads.** `copy_assign` writes no assignment row for a family with one copy, so
`--psv-clusters` has no reads to cluster there. Two consequences follow, both found while scoring:
- **The 48 two-copy held-out families are untested.** After masking they have one copy, and none of their masked-copy reads
  is placed. They are not eligible, so R1 and R2 stand on the 12 families with >= 3 copies.
- **R3 is not measurable.** The single-copy control (300 random genes and 26 `high_het`, each a one-copy family) placed no
  read, so its 0 / 300 false `copy` calls are vacuous and are not counted as a pass.

Amendment 1 was chosen on the understanding that a one-copy family is clustered like any other. That understanding came from a
check of mine that confirmed the run finished, not that it placed reads, so it was wrong.

#### Verdict (2026-10-15)

**Reading 3: not claimed.**
- R1 fails: 4 of 12 eligible masked copies are recovered held out, against a bar of 0.70; development gives 15 of 25.
- R2 fails as well.
- R3 cannot be measured.

#### Post hoc (descriptive, not pre-registered; `bench/o3_realdata/posthoc.py`)

- **About half the misses are fragmentation.** The masked copy's placed reads often sit in several PURE unexplained clusters
  (purity >= 0.9), none of which holds 80% of them alone. Taken together, those pure clusters hold >= 80% of the reads in
  8 of 12 held-out copies (4 by the registered single-cluster rule) and 18 of 25 development copies (15).
- The remaining misses:
  - **unlinked:** the reads form a cluster the linkage test sets aside (MCL20, GWFAM4);
  - **impure:** the reads share an unexplained cluster with survivor reads (MCL34, MCL69, GWFAM123, GWFAM402, GWFAM407,
    GWFAM415);
  - **explained or unresolved:** MCL27;
  - **in no cluster:** GWFAM47;
  - **just under the bar:** GWFAM331, with combined recall 0.795.
- **Specificity on real reads is poor.** In the added-back arm every copy is present, yet unexplained clusters with
  est_div > delta (`copy`) remain in 17 of 60 held-out families and 16 of 53 development families.

  | | held out | development |
  |---|---|---|
  | unexplained clusters labelled `copy` / `allele_or_copy` | 57 / 81 | 27 / 183 |
  | `copy` clusters made mostly of survivor reads | 41 of 57 | 21 of 27 |
  | `copy` clusters made mostly of the masked copy's reads | 16 | 6 |

  The simulation flagged 0 of 152 intact families. So on real reads, a `copy` cluster in a family whose copies are all in
  the genome is common.
- **What a fix would need** (for a new prereg, not done here):
  - place a one-copy family's reads on its copy when `--psv-clusters` is on, which gives 2-copy losses and the single-copy
    control a test;
  - decide whether fragments of one source should count together (they did in 4 held-out families);
  - find what the `copy` clusters of intact families are before relying on the label: collapsed copies in the assembly,
    isoforms or exons the catalog copy does not carry, or reads from other loci.

#### Follow-up (post hoc, 2026-10-15): what are the `copy` clusters of intact families?

The added-back arms have every copy in the genome, yet they hold 57 `copy`-labelled unexplained clusters (held out) and 27
(development). Each was traced (`bench/o3_realdata/copyclusters.py`, `copyclusters_summary.py`; outputs
`/mnt/linuxdisk/tmp/o3_realdata/{heldout2,dev}/classified.tsv`, `copyclusters_summary.json`) with four kinds of evidence:
- each read's copy of origin (its baseline primary);
- the prototype's clustering re-run on each host. It matched every Rust cluster exactly (Jaccard 1.0, 84 / 84), and splits
  each cluster's unexplained sites into those no catalog copy covers and those covered;
- the reads aligned to KB3781's own maternal and paternal assemblies, with the family's copies lifted onto both;
- for clusters matching neither haplotype, the reads aligned to human (CHM13), chimpanzee (mPanTro3) and orangutan
  (mPonPyg2).

Rules, applied in order:
1. **other family's copy**: most reads match a haplotype (>= 0.995 identity, >= 0.9 coverage) at another family's copy.
2. **not in the genome**: most reads match neither haplotype at >= 0.995. These are then given the species with the
   highest median identity if it is >= 0.99. Below that, they are gorilla-like when gorilla is closest by >= 0.02, else
   unresolved.
3. **missing from the reference**: the reads match a haplotype locus that is no copy of the family.
4. Own-family reads:
   - **catalog lacks the exons**: >= half the unexplained sites lie where no catalog copy has sequence;
   - **minority segment**: the sites are covered by a median < 0.2 of the cluster's reads;
   - **own-copy variant**: otherwise.

| what the cluster is | held out: clusters / reads | development: clusters / reads |
|---|---|---|
| reads of another gorilla family's copy | **35 / 2,292** | 0 / 0 |
| other apes' transcripts (chimpanzee; orangutan) | 4 / 103; 1 / 5 | 2 / 13; 2 / 8 |
| gorilla-like sequence in neither haplotype assembly | 3 / 1,985 | 0 / 0 |
| a copy missing from the reference assembly | 0 / 0 | **3 / 469** (GWFAM175) |
| the catalog copy lacks the exons | 4 / 320 | 2 / 16 |
| a few reads' divergent segment flips the cluster consensus | 0 / 0 | 5 / 930 |
| allele-like variant of the copy's own reads | 5 / 129 | 4 / 39 |
| unresolved (no species >= 0.99) | 5 / 183 | 9 / 38 |
| total | 57 / 5,017 | 27 / 1,513 |

**Reading.**
- **Reads from other loci** cause most of these clusters. On the shipped catalog, 35 of 57 are reads of another
  family's copy. They sit in 7 families (MCL69 10 clusters, MCL126 9, MCL43 6, MCL288 5, MCL296 3, MCL80 1, MCL286 1).
  Their reads match that other locus closely (each cluster's median read identity to its best haplotype hit: median
  0.999, range 0.994-1.0). `--no-as-tied-only` places them on this family, whose catalog lacks their source.
- **Other apes' transcripts** (chimpanzee, orangutan) remain in the library, as found before (A44).
- **Collapsed or missing copies are rare.**
  - Only GWFAM175's reads match a haplotype locus that is not a copy of the family (cluster medians 0.9988-0.9996): a real copy missing from
    the reference, and the case O3 is for.
  - MCL20's and MCL33's large clusters are gorilla sequence 2-3% from any locus in either haplotype. They may be copies in
    assembly gaps, but the assemblies cannot say.
- **Catalog and consensus artifacts** make the rest:
  - exons the catalog copy lacks are unexplained by construction;
  - a few reads with a locally divergent segment can turn a cluster's consensus to `alt` where only they cover the locus
    (GWFAM47: 6 of 176 reads, a ~70 bp window; the median read has 0 mismatches).

**What would cut the false `copy` calls** (each a method change, needing its own prereg):
- keep only reads whose primary or best alignment lies in the family, or cluster related families together;
- screen out other species (the A44 read-share screen);
- count only sites that a catalog copy covers;
- require the unexplained alleles to be carried by most of the cluster's reads.


---

## PREREG_o3_unexplained_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_unexplained_2026-10-15.md`)*

## PREREG — O3 (advisor's definition): PSV clustering of the reads O2 places on a family, held-out simulation (2026-10-15)

**Written after development on seeds 31, 32, 33 and before any held-out world is built.** Design:
`docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md`. Frozen code (md5): `bench/o3_unexplained/psvclust.py`
e99d49b6bd3d9d1ad9f61ae19fa70c98, `score.py` 852aa244a00967c153abb6f2b9bc2cc1, `world.py` dfe160287b5fe50e07a81487d48eb237;
`copy_assign` built `--release` from 5022350c (sha1 87eed3f8); minimap2 2.30, samtools 1.22.1.

### Claim under test

Given the reads of a multi-copy family mapped to a reference that LACKS one of its expressed sequences (a copy, or an allele), O2
places those reads on the remaining copies; clustering each copy's reads by read-derived variant sites (psvclust.py) recovers the
missing sequence's reads as one cluster that no reference copy explains, does not produce such a cluster when nothing is missing,
separates two missing sequences absorbed onto one copy, and labels the cluster `copy` when its estimated divergence exceeds
delta = 0.00958 (else `allele_or_copy`: from RNA alone an allele and a copy are the same signal below delta).

### Method (frozen; psvclust.py docstring)

Per host class (copies with an identical spliced sequence are one class), the reads whose O2 best copy (`catalog_copy_idx`) is in
it and every family copy's spliced sequence are aligned to the class's padded locus; sites = positions where >= 3 reads carry the
same non-reference base; reads join the compatible cluster with the most co-covered sites; clusters < 3 reads dissolve; clusters
merge unless a site separates them with >= 3 reads on each side; a cluster is UNEXPLAINED when >= 2 of its alt sites are carried by
no family copy (and no copy explains all of them); est_div = unexplained alt sites / bases covered.

### Held-out worlds and arms

`world.py build OUT --seed S --strands` for **S = 41, 42, 43, 44, 45**; read variants **E0** and **E2**; `world.py run` with arm
**P** (`copy_assign --families --flag-missing-copies --no-as-tied-only`: every molecule placed on the family) — the primary arm — and
arm **T** (the shipped exact-tie gate) reported beside; then `psvclust.py` and `score.py rows`. 51 families per world.

### Bars (arm P, pooled over 5 seeds × 2 variants)

| id | measure | bar | dev (s32-33, arm P) |
|---|---|---|---|
| **D1** | absent copies with d >= 0.6% recovered (pure unexplained cluster: >= 80% of the reads O2 placed, purity >= 0.9) | >= 0.90 | 112/112 |
| **D2** | control families (NULL_TWIN, NULL_MULTI, NULL_CLOSE, ADDBACK, DECOY) with any unexplained cluster | <= 2 of 120 | 0/48 |
| **D3** | COMBO_EH and DOUBLE families with BOTH missing sequences recovered as separate clusters | >= 0.80 | 8/8 |
| **D4** | recovered alleles labelled `allele_or_copy`; recovered copies with d >= 2% labelled `copy` | >= 0.90 each | 28/28; 60/60 |

### Reported beside (no bar)

Copies at 0.3% (below the origin certificate's absorption floor) and at 1% (at delta: either label); every scenario; arm T on every
measure; the existing detector (`o3_flag == missing_copy` at the host copy) on D1 and D2 as the baseline; est_div against d.

### Reading, fixed now

- **D1-D4 pass** → PSV clustering of the reads O2 places on a family is the O3 method in its advisor's definition, at simulation
  level: it finds missing copies down to the allele range, does not flag intact families, separates two missing sequences, and
  states which ones only DNA can tell from an allele. Next: the same with real families on the main machine — remove real copies from
  the GENOME, re-map the real reads, add them back — and a single-copy-gene control for the real allele background (pre-registered
  separately). The method is then ported to Rust inside `copy_assign` under its own spec.
- **D1 or D2 fails** → the method is not adopted; the failing scenarios are named; nothing is re-tuned on the held-out worlds.
- **D3 or D4 fails alone** → detection is claimed, separation / labelling is not.
- Arm T is not a bar: on dev it sees only the exact-tie (twin) scenarios (16/16 there, 0 elsewhere), which is itself the finding to
  report — the absorbed reads of a missing copy are mostly UNIQUE mappers, not ties.

### Known limits, stated before the run

Simulated reads carry uniform HiFi errors, equal expression and one or two absent sequences per family; real data add gene
conversion, many-copy families, real error structure, and a heterozygous allele background at every expressed gene of a diploid
individual, which this method reports as `allele_or_copy` clusters — the real-data specificity of the `copy` label is therefore the
open question for the next level.

## OUTCOME (2026-10-15) — D1-D4 all PASS on the held-out worlds

`pixi run bash bench/o3_unexplained/heldout.sh ~/rustle_data/o3_unexplained/heldout` (seeds 41-45 `--strands`, E0 / E2, arms P and
T); frozen code unchanged (md5s as above). Arm P:

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| D1 | absent copies >= 0.6% recovered as a pure unexplained cluster | **310 / 310 = 1.000** (existing detector 292/310) | >= 0.90 | ✅ PASS |
| D2 | control families with an unexplained cluster | **0 / 120** (existing detector 0/120) | <= 2 | ✅ PASS |
| D3 | two missing sequences on one host recovered separately | **40 / 40** | >= 0.80 | ✅ PASS |
| D4 | alleles labelled `allele_or_copy`; copies >= 2% labelled `copy` | **74 / 78 = 0.949; 170 / 170** | >= 0.90 each | ✅ PASS |

**Reported beside.**
- Below the origin certificate's absorption floor: copies at 0.3% recovered 40/40 (existing detector 4/40); alleles at 0.25% 18/20
  (the two misses, both seed 43, are the only unrecovered sources of the run).
- Labels at delta: copies at 1% split between `copy` and `allele_or_copy` (TWIN 10/10, UNIQ 14/6, MULTI 12/8, DOUBLE 6/14); the 4
  alleles labelled `copy` are all at 0.8%, next to delta. The label is a divergence estimate and is honest only away from delta.
- Low expression: copies at 8, 12 and 16 reads 30/30.
- **Arm T** (the shipped exact-tie gate, the advisor's literal "multi-mappers") recovers only the twin scenario (100/100 there, 0
  elsewhere; D1 80/310): the reads of a missing copy are exact ties only when the reference holds identical twins of it; otherwise they
  map UNIQUELY to the nearest copy with extra mismatches. O3 must therefore cluster every read O2 places on a family, not only ties.

**Reading (as fixed).** At simulation level, PSV clustering of the reads O2 places on a family is the O3 method in its advisor's
definition: it finds every missing copy from 0.3% upward, flags no intact family, separates two missing sequences on one copy, and
labels by divergence which ones RNA cannot tell from an allele. Next: real families on the main machine (remove real copies from the
GENOME, re-map the real reads, add them back) with a single-copy-gene control for the real allele background, pre-registered before
it runs; then the Rust port inside `copy_assign`.


---

## PREREG_o3_graph_round10_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round10_2026-10-16.md`)*

## Round 10: the integrated pipeline v3 under calibrated v2 reads, end to end, on an eleventh set of full families (2026-10-16)

Pre-registration, written before any world of the eleventh set was built: no read of these families has been simulated. Code:
branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal (2026-10-16): "integrate O1, O2 and O3 more so it's faster and cheaper to test full families; also
improve the simulations to account for more realistic ones and determine if PSVs can also include SUNs and probably some
influence of the cigar string, nm, de tags". Round 9 failed under both profiles. The development work since then removes
the round-9 failure mechanisms; see "Known before running" below and `docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 11.
The round-9 failures were:
- O2 left realistic reads `unexplained` on catalogued copies (1-3.5% fewer correct reads);
- O1 dropped expressed copies;
- O3 made false calls on present copies' reads.

**Question.** On a fresh set of full (uncapped) families, with reads from both calibrated v2 profiles, does the integrated
pipeline v3 do at least as well as the current method in every world type? That means at least as many correct read outcomes
and no more wrong ones.

### The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v3** (every change since round 9 was chosen on development sets 8-10, all spent):
- **O1:** `o1seed.py`. A curated copy is kept with >= 5 primaries OR >= 1 owned primary at MAPQ >= 1 (e40367ea).
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact`, then `o2graph.py run W integl 2`.
  - `--union-exact` (36a430f8): a paralog isoform projects onto a copy even when the spliced aligner drops a terminal exon
    < 30 bp. A projected structure stays apart from the copy's own unless its junctions match exactly.
  - `integl` (eb4f0c95), the genome fallback for the reads the paths leave `unexplained`:
    - every primary and secondary genome alignment is scored as mismatches + indel events;
    - an insertion or clip >= 30 bp refuses (it is O3's evidence);
    - the read's best alignment off the family's seeded copies keeps it `unexplained`;
    - the read is placed only at margin 2 within delta.
  - O3 still reads `o2.integ.m2.tsv`, so the fallback changes no O3 evidence.
  - The set is scored on `o2.integl.m2.tsv`.
- **O3:** `detect.py W --member-genome-check --member-attribution --arms I`.
  - `--member-attribution` (12f8f2af): an absent-copy PSV call keeps only the reads its consensus explains better than every
    single catalogued copy, by O2's margin of 2.
  - The mosaic rule's nearest copy compares off-site positions with what the cluster's reads show, not with the host's base
    (f080a30f).
  - Speedups with identical calls: lazy isoform maps, faster PSV clustering (9829bc06).

Both systems run per world from one command. The O2 and O3 speedups (`--o2-procs 2`, `--o3-procs 2`) leave outputs
byte-identical.
```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout10`, committed in 70ebf2b1):
- CYP3A, FAM47, FAM9, DAZ: round-9 candidates past its first 8, never built into a world;
- HNRNPCL, PRR20, ZNF705, TAF11L, F8A, FAM236, CENPVL, FOXD4L, OR2T, KRTAP10, H2BC, FAM66, CXADRP, GOLGA2P: none was in any
  earlier list.

**Build rule:** the locus link of the rebuilt simulation sets.
- A gene's canonical transcript, spliced onto the other gene's locus, covers >= 0.8 of it in one direction, at identity >= 0.90.
- Span <= 60 kb.
- **No copy cap** (`--cap 1000`).
- `families.py build families10_all.json --cap 1000 --link locus --min-cover 0.8`, every group.

The set is the first 8 families with >= 3 copies, in list order. CYP3A links 2 copies; FAM47 and FAM9 form no family.

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| DAZ | 3 | 7 | chrY |
| HNRNPCL | 3 | 3 | chr1 |
| PRR20 | 3 | 4 | chr13 |
| ZNF705 | 6 | 10 | chr8, chr11, chr12 |
| TAF11L | 45 | 45 | chr5 |
| F8A | 3 | 3 | chrX |
| FAM236 | 4 | 7 | chrX |
| CENPVL | 3 | 3 | chrX |

70 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families10.json` (this set only), `families10_all.json`, `families10.log`.

### Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under two read profiles:
`bench/o3_graph/profiles/A119b_testis_v2.json` and `bench/o3_graph/profiles/public_testis_v2.json`. That makes 12 worlds.
- **The v2 profiles (round-9 review):**
  - errors and donor alleles are taken from every MAPQ-60 read;
  - the expression spread is measured over read copies, less Poisson noise;
  - truncation is redrawn below 200 bases;
  - plus round 9's calibrated truncation, homopolymer-biased indels, het and hom donor alleles, and retained introns.
- With a profile, `--het 0` means the donor's alleles come from it.
- **The seed:** `precheck.py families10.json heldout10 121 131` (`precheck10.txt`) gives the smallest seed in 121-130 with
  >= 3 detectable removed copies: **seed 122** (42.6 s, 14.4 GB peak):
  - detectable: HNRNPCL 0.0453, TAF11L 0.0444, ZNF705 0.0244;
  - within delta: DAZ 0.0, PRR20 0.0011, F8A 0.0, FAM236 0.0022, CENPVL 0.0006.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families10.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout10
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout10 122 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout10 122 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```
`fastset.py` writes `e2e.json` (per-read detectability, on `o2.integl.m2.tsv`) and `bars.txt` (`summary.py bars`) per profile.
Full families with a 45-copy array may need copy_assign under the heavy lock. If a run's current step exceeds the light
limit, the run is repeated with `--heavy-current`, whose output does not depend on it.

### Read outcomes

These are round 9's (`summary.py` `read_fate`, `cmd_e2e`), unchanged:
- per-read detectability: a REMOVED copy's read is detectable by its own primary `de` on the masked genome (> delta);
- a BRANCH read is detectable by construction;
- an O3 call on a read takes precedence over its O2 placement;
- a read called in several families is judged by its call in its own family.

### Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither (`summary.py bars`).

### Reading

- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "with reads calibrated on two real libraries, the integrated O1 -> O2 ->
  O3 pipeline v3 gives, in every world type of a fresh set of full families, at least as many correct read outcomes as the
  current method and no more wrong ones".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.

### Known before running (development sets 8, 9 and 10, v2 profiles, per-read scoring)

Every development run below was simulated fresh with the v2 profiles and the new O1 rule:
- `fast/fbdev/` (O3 runs: `calls.tsv.prev`, `calls.noattr.tsv`, `calls.tsv`);
- `fast/fbdev_exact/` (`--union-exact`).

| configuration | set 8 A119b | set 8 public | set 9 A119b | set 9 public | set 10 A119b | set 10 public |
|---|---|---|---|---|---|---|
| no fallback (round 9's O2), old O3 | FAIL | FAIL | FAIL | PASS | FAIL | FAIL |
| fallback, old O3 | FAIL | FAIL | FAIL | PASS | FAIL | PASS |
| fallback, O3 consensus fix | FAIL | FAIL | FAIL | PASS | FAIL | PASS |
| fallback, + member attribution | PASS | FAIL | PASS | PASS | FAIL | PASS |
| **v3: fallback, attribution, `--union-exact`** | **PASS** | **PASS** | **PASS** | **PASS** | **FAIL** | **PASS** |

- **v3 fails only set 10 A119b, by 2 correct reads in INTACT** (1,136 vs 1,138; errors 0 vs 2). The lost reads are short
  reads whose sequencing edits exceed delta's share of their length: 2-5 edits on 184-484 bp, MAPQ 60, on their own copy.
  The current method assigns them.
- **The rejected alternatives,** measured on the same runs:
  - `combined` mode, the genome on every read with the elsewhere guard: errors drop (OUTSIDE 76 -> 0), but correct reads
    collapse in arrays (set 9 INTACT 935 -> 534).
  - `seeded` mode, the genome on explained reads without the guard: fixes CROSS but costs INTACT 13-47 reads.
  - `resolve` mode, the genome on O2's tied reads: adds 0-1 read per world.
- **The development summary of v3's margins** (correct, integrated minus current):
  - INTACT: +24, +97, +28, +45, -2, +6;
  - CROSS: +157, +230, +23, +46, +8, +13;
  - SILENT: +57, +176, +34, +41, +13, +11.
- **Errors:**
  - in INTACT, CROSS and SILENT: 0-1 per world, against 0-9 for the current method;
  - in REMOVED: 0-2 (current 0-12);
  - in OUTSIDE and BRANCH: far fewer than the current method.
- **A per-world bar is decided by single reads** wherever the current method makes no error, and by a handful of reads in
  INTACT, where both systems see the full catalog.

### Reported beside the bars (not gating)

- the fates per read class and the error sources;
- O1's drops and new copies;
- arm I's calls, true or false;
- copy-level scoring (`summary.py e2e` without `--per-read-detectability`);
- **the run's cost:** `timing.json`, every step's wall time and peak RSS;
- **an untrimmed-library run:** the public profile with adapters (`public_testis_v3.json`) on the same seed, trimmed by
  `fastset.py --trim-adapters profiles/public_testis_adapters.json`, same options, scored the same way.

### After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1308.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-16; run with the code of 58bca61d, started 19:30)

### Verdict: FAIL under both profiles. Not claimed.

Bars from `fastset.py`'s `e2e.json` (per-read detectability, on `o2.integl.m2.tsv`):

| world | A119b: correct (current / integrated) | errors | B1, B2 | public: correct | errors | B1, B2 |
|---|---|---|---|---|---|---|
| INTACT | 645 / 736 | 0 / **8** | yes, **no** | 576 / 654 | 0 / 0 | yes, yes |
| CROSS | 689 / 778 | 0 / **8** | yes, **no** | 616 / 694 | 0 / 0 | yes, yes |
| OUTSIDE | 828 / 830 | 142 / 61 | yes, yes | 802 / **751** | 115 / 57 | **no**, yes |
| REMOVED | 879 / 959 | 8 / 8 | yes, yes | 843 / 919 | 6 / 0 | yes, yes |
| BRANCH | 682 / 1,043 | 154 / 65 | yes, yes | 576 / 946 | 139 / 48 | yes, yes |
| SILENT | 464 / 475 | 0 / **8** | yes, **no** | 424 / 439 | 0 / 0 | yes, yes |

B1 holds in 11 of 12 worlds, by +2 to +370 correct reads. The integrated pipeline makes fewer errors in OUTSIDE, BRANCH
and REMOVED. Copy-level scoring (`e2e_copylevel.json`) gives the same verdict; there it also fails REMOVED's B2 (A119b 9 vs
0, public 2 vs 1).

**The two failures, each one mechanism, found post hoc:**
1. **A119b, B2 in INTACT, CROSS and SILENT: one false O3 call, in TAF11L, the 45-copy array.**
   - Arm I calls 8 reads of TAF11L.33, a present copy, an absent copy, on host 44 (est_div 0.012-0.016, members' median
     genome `de` 0.0103).
   - The same call recurs in BRANCH, where it costs nothing.
   - **Cause:** the family's unexplained pool (281 reads from about 40 copies) clusters on one host into a few large mixed
     clusters. Their consensus carries several copies' private alleles, so it is `unexplained` beyond delta.
   - Member attribution then keeps the TAF11L.33 reads whose own alleles match that consensus at >= 2 positions where
     TAF11L.33's catalogued sequence differs.
   - Development's arrays (set 9's USP17L and TSPY) did not produce such a call.
2. **public, B1 in OUTSIDE (-51): O1's new copy within delta of a catalogued one.**
   - O1 seeds the uncatalogued FAM236.1 correctly, and the integrated pipeline gets 40 reads of uncatalogued copies right.
   - FAM236.1 differs from FAM236.0 by about 0.2% (the pre-check's 0.0022), so 166 catalogued FAM236 reads become `thin` or
     `tied` between the two.
   - The current method, which does not know FAM236.1, assigns those reads, and makes 114 errors on uncatalogued reads.

**Reported beside the bars:**
- **Arm I:**
  - REMOVED: all 3 detectable removed copies called under both profiles (HNRNPCL.1, TAF11L.41, ZNF705.0);
  - BRANCH: 5 of 5 novel-exon copies under both profiles;
  - false calls: the TAF11L.33 call above (A119b: INTACT, CROSS, BRANCH, SILENT); none under the public profile.
- **O1:** 5 silent copies dropped per SILENT world and no expressed copy; 3-4 new copies in OUTSIDE.
- **Cost:**
  - A119b 327 s for six worlds, including the one-time local reference build (44 s, 14.4 GB);
  - public 273 s, peak 1.9 GB (copy_assign);
  - per set: O3 170-185 s, current 131-144 s, O1 53-61 s, O2 33-40 s, the genome fallback 2-3 s.
- **The untrimmed-library run:** below.

**The untrimmed-library run** (secondary; the public profile with its Tn5/SMART adapters, `public_testis_v3.json`, trimmed by
`--trim-adapters`, same seed and options): FAIL. Correct reads match the adapter-free public run (INTACT 656 vs 654, BRANCH
947 vs 946), so trimming restores the reads; two known mechanisms add errors:
- **INTACT, CROSS, SILENT, 2 errors each:** O1 dropped PRR20.1, an expressed copy all of whose reads tie, and O2 placed two of
  its reads on PRR20.0.
- **REMOVED, 17 errors:** a second mixed-cluster call in TAF11L (host 43, 25 members from present copies), round 10's
  failure 1.
- **OUTSIDE B1:** 753 vs 799, failure 2.
- Its O3 step ran while `detect.py` gained the opt-in `--member-attribution-delta`. The run did not pass it, so arm I behaved as
  at 58bca61d.

### Reading

**The pre-registered claim fails under both profiles. No partial claim is made.** Set beside round 9, on held-out families:
- **The catalogued-read loss under realistic reads is gone.** B1 holds in INTACT, CROSS and SILENT under both profiles, by +11 to
  +91 reads; round 9 lost 1-3.5% there.
- **The remaining errors are O3 calls in a 45-copy near-identical array** (A119b), and one O1 new copy within delta of a
  catalogued copy that makes that copy's reads ambiguous (public OUTSIDE).
- Neither mechanism occurred on the development sets.

### Independent review (2026-10-16, 19:53-20:21) and corrections

Fresh subagent; full report `/mnt/linuxdisk/tmp/o3_graph/review10/REVIEW.md`.
- **The verdict holds.** Its own scorer, written from this document, reproduces `e2e.json`, `bars.txt`,
  `bars_copylevel.txt` and the untrimmed run's bars exactly.
- **58bca61d reproduces byte-identically:** O1, O2 and the fallback in three worlds, and arm I in three.
- **The protocol is clean:**
  - the commit precedes every world;
  - none of the 344 earlier worlds holds a read of these families;
  - the selection and the seed follow their rules;
  - the commands are the pre-registered ones, and all 130 steps exited 0.

**Corrections to the results above:**
1. **Failure 1's mechanism was wrong (Critical).** The facts stand: one arm-I call, 8 TAF11L.33 reads, host 44. The cause is
   different.
   - **The host window holds a neighbour's fragment.** TAF11L.44's +-2 kb window contains the last 225 bp of TAF11L.43,
     1,775 bp away. It is the only pair of catalogued copies in the set closer than the pad; no development host had one,
     not even the full 59-copy USP17L.
   - **The reads aligned to it are compared with nothing.** 13 haplotype-2 TAF11L.33 reads align to that fragment, which
     no catalogued copy's map covers. Their differences there count as unexplained alts; 8 of the 9 are TAF11L.33's own
     reference bases.
   - **Attribution kept them by default.** It keeps a read when no copy covers half of its positions.
   - The genome check passed the call by 0.0007 (median `de` 0.0103 over the 13 attribution-kept reads; the per-read check
     then kept 8).
   - The consensus did not carry other copies' private alleles. 038890fd (`--member-attribution-delta`), built on that
     wrong mechanism, does not remove the call.
   - **The fix is 88c445d9:** attribution keeps no read that no catalogued copy competes for.
2. **The same 8 reads are A119b REMOVED's 8 errors (Important).** There they sit inside the true TAF11L.41 call on host 43.
   REMOVED is a tie (8 vs 8), not "fewer errors".
3. **Failure 2's twin is FAM236.2, not FAM236.0 (Important).** It is 1 mismatch over 451 bp, `de` 0.0022.
   - The 166 lost reads are FAM236.2's.
   - Seeding FAM236.1 won 0 of its own reads (public) and 1 (A119b). The 40 uncatalogued correct reads are HNRNPCL.1's (24)
     and PRR20.2's (16).
   - The same mechanism costs 152 reads in A119b OUTSIDE, which passes B1 by only +2.
4. **Not disclosed (Important):** arm I takes its read hosts and hosts from copy_assign's run (`cur/A.assignments.tsv`). So
   the integrated pipeline still needs the current method's run (131-144 s per set, about 2 GB).
5. **Smaller points:**
   - **The untrimmed TAF11L call:** 17 present-copy reads and 8 correct TAF11L.41 reads.
   - **The PRR20.1 drop** also errs in OUTSIDE and BRANCH, through a MAPQ-0 tie-break.
   - **OUTSIDE's gross catalogued loss** is 161 and 169 reads, so "the catalogued-read loss is gone" holds for INTACT,
     CROSS and SILENT only.
   - **Remaining errors also include** OUTSIDE's CENPVL.1 reads (61 and 57) and BRANCH reads lacking the novel exon (57 and
     48).
   - **O1 kept 18 of 23 silent copies.**
   - **"Development set 10's 2 lost reads" is a net figure:** 19 lost and 17 gained. 4 of the 19 are `thin` at MAPQ 10-30.
   - **The OUTSIDE B1 count depends on the scoring rule:** counting O3 flags of uncatalogued in-genome copies as correct
     would fail A119b OUTSIDE B1.
   - **The local reference** is still unvalidated for v2 reads (round 9's caveat).


---

## PREREG_o3_graph_round11_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round11_2026-10-16.md`)*

## Round 11: the integrated pipeline v4 under calibrated v2 reads, end to end, on a twelfth set of full families (2026-10-16)

Pre-registration, written before any world of the twelfth set was built: no read of these families has been simulated. Code:
branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal is unchanged: integrate O1, O2 and O3 so full families are fast and cheap to test, make the
simulations more realistic, and determine what PSVs/SUNs and the CIGAR, NM and `de` tags contribute.

Round 10 (register 1308) failed. Its independent review corrected the first failure's mechanism to a host-window artefact. v4
removes each round-10 failure mechanism found on development data (`docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 12).

**Question.** On a fresh set of full (uncapped) families, with reads from both calibrated v2 profiles, does the integrated
pipeline v4 do at least as well as the current method in every world type? That means at least as many correct read outcomes
and no more wrong ones.

### The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v4** = round 10's v3 plus three changes, all chosen on development data (sets 8-10 and round 10's spent worlds):
- **O1:** `o1seed.py W --twins`.
  - **The keep rule:** >= 5 primaries, or >= 1 owned MAPQ >= 1 primary.
  - **A new copy within delta of a kept curated copy is its twin** (genomic spans, minimap2 asm20 `de`).
  - **A curated copy the rule would drop is kept when it lies within delta of a kept one.** The reads cannot say which of the
    two is expressed.
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `run W integl 2 --error-floor`.
  - This is round 10's union paths and genome fallback.
  - **A twin takes a read from its curated copy only by the margin.**
  - **The `unexplained` gate is delta plus the world's error floor.**
- **O3:** `detect.py W --member-genome-check --member-attribution --error-floor --arms I`.
  - **Attribution keeps no read that no catalogued copy competes for.**
  - **The genome check subtracts the error floor.**
  - **The read hosts still come from copy_assign's run** (`cur/A.assignments.tsv`), as in rounds 7-10. O1 hosting
    (`--o1-hosts`) was measured and not adopted: it loses REMOVED and BRANCH recall in arrays.
- **The error floor** is the world's median `de` over MAPQ-60 primary alignments, i.e. its reads' sequencing error. Delta
  (0.00958) is the p99 exonic divergence between KB3781's haplotype assemblies (register 1199), an error-free divergence.
  Every read-level comparison with delta therefore subtracts the floor first.

Both systems run per world from one command:
```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --twins --error-floor --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```
`--o2-procs 2` and `--o3-procs 2` (defaults) leave outputs byte-identical.

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout11`, committed in c034cc5d):
- FOXD4L, OR2T, KRTAP10, H2BC, FAM66, CXADRP, GOLGA2P: round-10 candidates past its first 8, never built into a world;
- H3C, H2AC, KRTAP2, KRTAP1, PRR23, ETD, SULT1C, UGT2A, CDRT15, NAIP, HLA-DRB, TCEAL, ARMCX, OR4C, OR4K, OR11H, SEC14L, DUX4:
  none was in any earlier list.

**Build rule:** round 10's (locus link >= 0.8 at identity >= 0.90, span <= 60 kb, no copy cap):
`families.py build families11_all.json --cap 1000 --link locus --min-cover 0.8`, every group. The set is the first 8 families
with >= 3 copies, in list order (GOLGA2P links 2 copies):

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| FOXD4L | 6 | 6 | chr2, chr9 |
| OR2T | 3 | 3 | chr1 |
| KRTAP10 | 9 | 10 | chr21 |
| H2BC | 6 | 6 | chr6 |
| FAM66 | 3 | 10 | chr8, chr12 |
| CXADRP | 3 | 3 | chr15, chr18, chr21 |
| H3C | 3 | 3 | chr1 |
| H2AC | 4 | 4 | chr6 |

37 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families11.json` (this set only), `families11_all.json`, `families11.log`.

### Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under round 10's two read
profiles, `profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json` (`--het 0`). That makes 12 worlds.
- **The seed:** `precheck.py families11.json heldout11 141 151` (`precheck11.txt`; 93.9 s, 13.2 GB) gives the smallest seed
  in 141-150 with >= 3 detectable removed copies: **seed 141**.
  - Detectable: H2AC 0.1027, H2BC 0.0943, KRTAP10 0.0544, FAM66 0.0448, FOXD4L 0.0153, CXADRP 0.0113.
  - Within delta: OR2T 0.0061, H3C 0.0.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families11.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout11
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --twins --error-floor --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout11 141 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout11 141 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```
`fastset.py` writes `e2e.json` and `bars.txt` per profile. The scoring is per-read detectability with the error floor, on
`o2.integl.m2.tsv`. If a run's current step exceeds the light limit, the run is repeated with `--heavy-current`, whose output
does not depend on it.

### Read outcomes

**These are round 10's, with one change:** a REMOVED copy's read is detectable when its own primary `de`, minus its world's
error floor, exceeds delta (`summary.py e2e --error-floor`).
- **Why:**
  - Delta is an error-free haplotype divergence.
  - A read's `de` adds its own sequencing errors.
  - Without the correction, a read whose copy lies within delta of a present copy counts as detectable merely through its
    errors.
- **The change applies to both systems.** It moves REMOVED reads with `de` in (delta, delta + floor] from detectable to
  undetectable: placing them becomes `unknowable` instead of an error, and an O3 call on them stays correct.
- **The rule was decided on development data after round 10.**
- **Secondary, not gating:** round 10's rule (`e2e_r10rule.json`, same command without `--error-floor`) and copy-level
  scoring are both reported.

Unchanged:
- a BRANCH read is detectable by construction;
- an O3 call on a read takes precedence over its O2 placement;
- a read called in several families is judged by its call in its own family.

### Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither.

### Reading

- **PASS** (B1 and B2 hold in all 12 worlds under the primary scoring). Claim: "with reads calibrated on two real libraries,
  the integrated O1 -> O2 -> O3 pipeline v4 gives, in every world type of a fresh set of full families, at least as many
  correct read outcomes as the current method and no more wrong ones, with per-read detectability corrected for sequencing
  error".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.
- **If the primary passes and round 10's rule fails,** the claim is stated with that qualifier.

### Known before running (development, v2 profiles)

Nine runs: development sets 8, 9 and 10 under both profiles, plus round 10's two profiles and its untrimmed-library run
(`fast/r11dev6/`):
- **Primary scoring:** v4 PASSES all 9.
- **Round 10's rule:** v4 passes 5 of 9. The 4 failures are all REMOVED B2: set 8 public 9 vs 7, set 9 A119b 1 vs 0, set 10
  A119b 33 vs 10, set 10 public 13 vs 4. Placements of reads with `de` between delta and delta + floor.
- **v4's margins (correct, integrated minus current):**
  - INTACT: +92, +78, +54, +101, +37, +48, +8, +8;
  - CROSS: +90, +78, +191, +234, +36, +49, +17, +15;
  - SILENT: +12, +15, +69, +114, +32, +42, +11, +7;
  - REMOVED: +90, +76, +67, +87, +20, +12, +1, +55.
- **Errors in INTACT, CROSS and SILENT:** 0-1 per world, against 0-9 for the current method.
- **The thinnest margins:** set 10 A119b INTACT +8 and REMOVED +1, against 1 error vs 4.
- **Untested in development:**
  - a family set without arrays or near-identical copies;
  - histone clusters (H2BC, H3C, H2AC), whose copies are short and single-exon.

### Reported beside the bars (not gating)

- the fates per read class and the error sources, tabulated;
- O1's drops (of how many silent copies), its kept-undecidable copies, its new copies and twins;
- arm I's calls, true or false, in every world;
- round 10's scoring and copy-level scoring;
- **the run's cost:** `timing.json`, every step's wall time and peak RSS; the integrated pipeline's own steps (O1, O2,
  O2_genome, O3) apart from the current method's;
- **an untrimmed-library run:** `public_testis_v3.json` + `--trim-adapters profiles/public_testis_adapters.json`, same seed
  and options.

### After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1309.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-16; run with the code of fe73bc8c, started 20:57)

### Verdict: PASS under both profiles. Claimed.

Bars from `fastset.py`'s `e2e.json` (primary scoring: per-read detectability less the error floor, on `o2.integl.m2.tsv`):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 1,728 / 1,917 | 186 / 0 | 1,779 / 1,984 | 159 / 0 |
| CROSS | 1,793 / 1,984 | 189 / 0 | 1,834 / 2,042 | 161 / 0 |
| OUTSIDE | 1,415 / 1,813 | 367 / 44 | 1,503 / 1,893 | 297 / 37 |
| REMOVED | 1,545 / 1,844 | 198 / 0 | 1,626 / 1,928 | 162 / 1 |
| BRANCH | 1,758 / 2,040 | 313 / 118 | 1,791 / 2,107 | 297 / 119 |
| SILENT | 700 / 857 | 154 / 0 | 603 / 753 | 138 / 0 |

B1 and B2 hold in all 12 worlds. **Round 10's scoring and copy-level scoring give the same verdict, PASS in all 12 worlds**
(`bars_r10rule.txt`, `bars_copylevel.txt`), so the claim does not rest on the scoring change.

**What the margins are made of (post hoc, and it qualifies the margins, not the verdict):**
- **Most of the current method's errors in INTACT, CROSS and SILENT are its own O3 calls** in the interleaved chr6 histone
  clusters.
  - H2AC and H2BC genes alternate there. copy_assign's PSV clusters mix the two families' reads (for example, an H2AC host
    cluster of 22 H2BC.2 and 2 H2AC.0 reads).
  - In A119b INTACT, 8 such calls hold 229 present-copy reads.
  - The integrated pipeline makes no call there.
- **Rescored without H2AC and H2BC** (both systems, both profiles), the verdict is still PASS in all 12 worlds, with smaller
  margins (`/mnt/linuxdisk/tmp/o3_graph/fast/r11_nohist/`):

| world | A119b correct (cur / int) | errors | public correct | errors |
|---|---|---|---|---|
| INTACT | 1,600 / 1,608 | 0 / 0 | 1,672 / 1,724 | 4 / 0 |
| CROSS | 1,665 / 1,675 | 3 / 0 | 1,727 / 1,782 | 6 / 0 |
| OUTSIDE | 1,304 / 1,551 | 182 / 44 | 1,395 / 1,669 | 158 / 37 |
| REMOVED | 1,396 / 1,542 | 47 / 0 | 1,493 / 1,671 | 44 / 1 |
| BRANCH | 1,614 / 1,717 | 120 / 109 | 1,673 / 1,834 | 135 / 114 |
| SILENT | 645 / 652 | 0 / 0 | 564 / 582 | 4 / 0 |

**Fates per read class** (`e2e.json`; abst = abstain, corr = correct, erro = error, unkn = unknowable):

| world | system | catalogued (A119b) | uncatalogued | absent, detectable | absent, undetectable |
|---|---|---|---|---|---|
| INTACT | current | abst 132, corr 1,728, erro 186 | | | |
| INTACT | integrated | abst 129, corr 1,917 | | | |
| OUTSIDE | current | abst 37, corr 1,415, erro 159 | abst 227, erro 208 | | |
| OUTSIDE | integrated | abst 49, corr 1,562 | abst 140, corr 251, erro 44 | | |
| REMOVED | current | abst 37, corr 1,399, erro 175 | | abst 147, corr 144, erro 23 | abst 15, corr 2, unkn 104 |
| REMOVED | integrated | abst 44, corr 1,567 | | abst 37, corr 277 | abst 7, unkn 114 |
| BRANCH | current | abst 132, corr 1,728, erro 186 | | abst 103, corr 30, erro 127 | |
| BRANCH | integrated | abst 129, corr 1,917 | | abst 19, corr 123, erro 118 | |

The public profile has the same shape (`public/e2e.json`). The integrated pipeline's remaining errors are:
- uncatalogued reads placed on a paralog (OUTSIDE 44 / 37);
- BRANCH reads lacking the novel exon, which count as detectable by construction (118 / 119);
- one REMOVED read (public).

**Reported beside the bars:**
- **Arm I:**
  - **no call in any world of either profile contains a present-copy read;**
  - **REMOVED:** calls for all 6 detectable removed copies (CXADRP, FAM66, FOXD4L, H2AC, H2BC, KRTAP10), plus a branch call
    for H3C's removed copy, which lies within delta; all members are the removed copies' reads;
  - **BRANCH:** 5 novel-exon copies (CXADRP, FAM66, H2BC, KRTAP10, OR2T).
    - A119b calls 4 of them: CXADRP, FAM66 and H2BC by branch calls, KRTAP10 by two 3-read PSV calls.
    - public calls 3: CXADRP, FAM66, H2BC.
    - OR2T's (21 and 15 reads) is missed under both.
- **O1:**
  - **SILENT:** drops 7 of 12 silent copies under both profiles, and keeps 1 or 2 FOXD4L copies as undecidable (within delta
    of a kept copy);
  - **no expressed copy dropped;**
  - **OUTSIDE:** 4 new copies (CXADRP, FAM66, KRTAP10, OR2T), none a twin.
- **Cost:**
  - **A119b:** 211.6 s for six worlds, including the one-time local reference build (14.4 GB);
  - **public:** 144.7 s, peak 1.58 GB;
  - **the integrated pipeline's own steps:** 89-90 s per set (O1 4 s, O2 6-9 s, the genome fallback 1 s, O3 77-79 s);
  - **the current method:** 55-64 s.
- **The untrimmed-library run** (public v3 adapters, trimmed): PASS in all 6 worlds (INTACT 1,804 / 1,983, errors 133 / 0;
  REMOVED 1,629 / 1,927, errors 164 / 1), 155 s.

### Independent review (2026-10-16, 21:08-22:10) and corrections

Fresh subagent; full report `/mnt/linuxdisk/tmp/o3_graph/review11/REVIEW.md`.

**The PASS holds.**
- The reviewer's own scorer reproduces `e2e.json`, `e2e_r10rule.json`, `e2e_copylevel.json`, the untrimmed run and the
  histone-free rescoring, cell for cell.
- **The variants also pass:** round 10's rule, copy level, O3 judged only in the read's own family, and without H2AC/H2BC.
- **The protocol is clean:**
  - the commit precedes every world;
  - none of 3,447 earlier truth files holds a read of these families;
  - the commands are the pre-registered ones;
  - no code changed until after the runs.
- **No leak.** O1, O2 and O3 read no truth file, and the error floor comes from the BAM.
- **Byte-identical reproduction from fe73bc8c:** O1, O2 and the fallback in 18 of 18 worlds; the current method in 7 of 7.

**Corrections and qualifications:**
1. **O3 with `--o3-procs 2` is not deterministic (Important).**
   - `detect.py` opens one FASTA handle before forking the family workers, which then share its file offset. A direct test
     returned about 100 wrong slices in 3,000 fetches.
   - A119b REMOVED's published arm-I output was not reproduced in 20 reruns. Integrated REMOVED correct varies 1,844-1,846
     (A119b) and 1,926-1,928 (public); errors never change, and no variant holds a present-copy read.
   - The verdict is unaffected, but the pre-registration's "byte-identical" is false for O3. The development claim of
     4efa6d33 missed the race. **Fixed before round 12** (one handle per process, as O2 already had).
2. **Arm I's REMOVED calls were misreported (Important).**
   - **The "branch call for H3C's removed copy" holds H2BC.2's reads.** H3C.1, identical to H3C.2, is never called, which is
     right.
   - **An H2AC-family branch call also holds H2BC.2's reads.** Those reads are also in H2BC-family calls, so the read-level
     bars do not move. At family level these are false calls.
3. **The histone account (Important).**
   - **All, not most,** of the current method's INTACT, CROSS and SILENT errors are its O3 calls: 96-100% in the histone
     clusters, the rest KRTAP10.
   - **The 8 calls hold 186 distinct reads, not 229 slots.**
   - **The mechanism** is round 10's host-window artefact, in the current method. Other-family reads enter a host's pool
     through MAPQ-0 secondaries over the window, and align there to UNCATALOGUED histone genes that no catalogued copy's map
     covers. The families hold 4 of 16 H2AC, 6 of 17 H2BC and 3 of 13 H3C genes.
   - **H3C is hit too.** Without H2AC, H2BC and H3C, the primary scoring and round 10's rule still pass. Copy level fails
     public REMOVED B2 (3 vs 7) and the untrimmed run's (3 vs 8): CXADRP.1 reads placed by the wider delta + floor gate. So the
     commit message and the overview's "holds without the histone clusters" are overstated; "without H2AC and H2BC" is
     accurate.
4. **The margins are the current method's O3 artefact, not an O2 gain (Important).**
   - With O3 removed from both systems (context only, not pre-registered), A119b is INTACT +3, CROSS +5, SILENT +3, REMOVED -5
     (errors 40 vs 0), and public +50, +52, +16, +32.
   - Under A119b the integrated O1 -> O2 is at parity with copy_assign's O2. The margins above that come from O3.
5. **The integrated pipeline's zero errors there rest on one filter (Important).** Before the genome check, arm I forms 12
   calls on present reads per INTACT world (the same histone cross-family clusters, and KRTAP10). Attribution removes none.
   The call-level genome check removes them all: median `de` at most 0.0073 against 0.0126 / 0.0113, a margin of about 0.005.
6. **OR2T's BRANCH miss is a pipeline failure (Important).** `vg msga` failed in every order for OR2T in every world that
   catalogues its three copies, so the branch arm was off for it.
7. **Minor points:**
   - **BRANCH errors:** in A119b, 117 of 118 lack the novel exon and 1 carries it.
   - **O1 seeds 4 of the 8 uncatalogued copies.** A single-exon uncatalogued copy can never be seeded (`observed_chains`
     skips unspliced reads): FOXD4L.0 and H2AC.0. H3C.1's reads, identical to H3C.2's, are all of the integrated OUTSIDE errors.
   - **Cost:** v4's arm I takes its read hosts from copy_assign's run, so the integrated pipeline costs about 145-154 s per
     set including that run, against 55-64 s for the current method.
   - **The fates table** covered only part of the worlds (full tables in the review).
   - **The scorer counts an O3 call in any family;** the verdict holds when O3 is judged only in the read's own family.
   - **Under the primary scoring, the REMOVED bar is close to tautological for the integrated pipeline:** its gate, its
     per-member check and the scorer share one threshold. The informative views are round 10's rule and copy level, and both
     pass.
   - **The local reference** is still unvalidated for v2 reads.


---

## PREREG_o3_graph_round12_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round12_2026-10-16.md`)*

## Round 12: the integrated pipeline v5, independent of the current method, on a thirteenth set of full families (2026-10-16)

Pre-registration, written before any world of the thirteenth set was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal is to integrate O1, O2 and O3 so that full families are fast and cheap to test. Round 11 (register
1309) passed with v4, whose O3 still took its read hosts from the current method's run (copy_assign,
`cur/A.assignments.tsv`). v5 hosts O3's reads from O1, so the integrated pipeline O1 -> O2 -> O3 needs nothing from the
current method; the current method runs only as the comparison.

**Question.** On a fresh set of full families, including an 80-copy array, does v5 do at least as well as the current method
in every world type, with both calibrated v2 profiles? That means at least as many correct read outcomes and no more wrong
ones.

### The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v5** = round 11's v4, with O3 hosted by O1 (`--o1-hosts`; f45e0fe5, f4f2b248):
- **O1:** `o1seed.py W --twins` (round 11).
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `run W integl 2 --error-floor`
  (round 11).
- **O3:** `detect.py W --member-genome-check --member-attribution --error-floor --o1-hosts --arms I`.
  - **Hosts:** O1's seeded copies, new ones included and dropped ones excluded.
  - **Reads:** every O1 family read is hosted on the seeded copy its primary alignment overlaps most. Failing that, on the
    copy a secondary alignment overlaps most; failing that, on the copy whose +-1 kb window its alignments overlap.
  - **No file of the current method's run is read.**
- **Since round 11's review:**
  - O3 opens one FASTA handle per process (7fbc236d; review I1). Two-worker and one-worker runs are byte-identical.
  - `--mapped-sites` (981708b5) was measured and is not adopted. It changes no bar and leaves the pre-check false calls as
    they were: these are cross-family branch calls and identical-copy clusters, all removed by the genome check.
- **Development: v5 passes all 12 runs under the primary scoring** (`fast/r12dev_base/`).
  - **The runs:** sets 8, 9 and 10 under both profiles, round 10's three runs and round 11's three, all spent.
  - **Round 10's rule:** v5 passes 8 of 12. The 4 failures are REMOVED B2 in sets 8 public, 9 A119b, 10 A119b and 10
    public, as for v4.
  - **Against v4 (copy_assign hosting), 72 world bars are identical in 68.** The four that differ:
    - round 11 A119b REMOVED 1,844 -> 1,845 (v4's run hit the FASTA race) and BRANCH 2,040 -> 2,041;
    - round 10 A119b BRANCH 1,043 -> 1,042;
    - set 8 A119b BRANCH 2,536 -> 2,538.
  - The -25 BRANCH reads seen earlier in set 9 public (`fast/r11dev7/`) do not recur with one FASTA handle per process.
  - Before the window fallback (f4f2b248), hosting only primaries lost array reads and passed 6 of 9.

The command:
```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --het 0 --jobs 2 --threads 2
    --sim-profile PROFILE
```
On development set 10 (A119b) this one command passes. The integrated steps take 118 s for six worlds (O1 7 s, O2 27 s,
genome fallback 1 s, O3 83 s), and the current method 47 s (`fast/r12dev_fastset/`).

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout12`, committed in 35bf1d95):
- KRTAP2, KRTAP1, ETD, CDRT15, NAIP, HLA-DRB, OR11H, DUX4: round-11 candidates past its first 8, never built into a world;
- POLR2J, KRTAP19, KRTAP3, MRGPRX: in no earlier list.

**Build rule:** round 11's, `families.py build families12_all.json --cap 1000 --link locus --min-cover 0.8`, every group. The
set is the first 8 families with >= 3 copies, in list order:

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| KRTAP2 | 4 | 5 | chr17 |
| KRTAP1 | 3 | 3 | chr17 |
| ETD | 3 | 5 | chrX |
| CDRT15 | 4 | 6 | chr17 |
| NAIP | 4 | 8 | chr5 |
| HLA-DRB | 3 | 3 | chr6 |
| OR11H | 3 | 3 | chr14, chr22 |
| DUX4 | 80 | 85 | chr1, chr3, chr4, chr9, chr10, chr16, chr20, chrY |

104 copies. DUX4's copies sit in D4Z4 tandem arrays, so neighbouring copies fall inside one another's +-2 kb host windows
(round 10's artefact). Files: `/mnt/linuxdisk/tmp/o3_graph/families12.json` (this set only), `families12_all.json`,
`families12.log`.

### Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under round 11's two
profiles, `profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json` (`--het 0`). That makes 12 worlds.

**The seed:** `precheck.py families12.json heldout12 161 171` (`precheck12.txt`; 76.5 s, 13.9 GB) gives the smallest seed in
161-170 with >= 3 detectable removed copies: **seed 161**.
- Detectable: KRTAP1 0.0514, HLA-DRB 0.0498, DUX4 0.0275, KRTAP2 0.0129.
- Within delta: CDRT15 0.0091, OR11H 0.0061, ETD 0.0032, NAIP 0.0.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families12.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout12
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --heavy-current --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout12 161 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout12 161 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```
- **Copy_assign's memory:** the current method runs under the heavy lock (`--heavy-current`) from the start. An 80-copy
  array may exceed the light jobs' 2 GB, and its output does not depend on the lock.
- **An 8-minute step timeout:** if a step exceeds rlock's 590 s timeout, the run is repeated with that step's timeout raised.
  The outputs do not depend on it. Either repeat is an amendment, recorded before any result is read.

### Read outcomes

Round 11's, unchanged:
- per-read detectability less the world's error floor (`summary.py e2e --error-floor`);
- a BRANCH read is detectable by construction;
- an O3 call on a read takes precedence over its O2 placement;
- a read called in several families is judged by its call in its own family.

**Secondary:** round 10's rule without the floor, and copy-level scoring.

### Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither.

### Reading

- **PASS** (B1 and B2 hold in all 12 worlds under the primary scoring). Claim: "the integrated O1 -> O2 -> O3 pipeline v5, which
  uses nothing from the current method, gives in every world type of a fresh set of full families, including an 80-copy
  array, at least as many correct read outcomes as the current method and no more wrong ones, with reads calibrated on two
  real libraries".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.

### Reported beside the bars (not gating)

- fates per read class, tabulated, and the error sources;
- O1's drops (of how many silent copies), kept-undecidable copies, new copies and twins;
- arm I's calls, true or false, in every world;
- round 10's scoring and copy-level scoring;
- **margins with each family left out in turn**, so that no single family carries the verdict (round 11's histone
  clusters did, for the margins);
- **the run's cost:** `timing.json`, with the integrated pipeline's own steps apart from the current method's.

### After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1310.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-16; run with the code of 00c241db, started 22:29)

### Amendment 1 (23:12, recorded after the A119b run's result was read)

**The public run crashed in its OUTSIDE world's O3 step.**
- The step outlasted rlock's 590 s limit (DUX4, 80 copies) and was killed before `/usr/bin/time` wrote its line.
- `fastset.py` then failed to parse the empty timing file and stopped (exit 1). No `e2e.json` or `bars.txt` was written for
  the public profile.

**As pre-registered for a step over the limit, the public run is repeated with the limit raised:**
- `--step-timeout 1800`, a new driver option; outputs do not depend on it;
- with the driver's parsing fixed (`parse_time`), so a killed step is reported, not crashed on;
- into `heldout12/public_r2`, with the crashed run's directory kept.

**This amendment comes after the A119b result was read,** which the pre-registration's amendment rule did not foresee.
A119b already FAILs (below), so the round's verdict is FAIL whatever the public repeat gives. The repeat completes the
record and is reported as such.

### Verdict: FAIL under both profiles. Not claimed.

Bars from `e2e.json` (primary scoring), A119b and the public repeat (`public_r2`, Amendment 1):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 2,189 / 2,380 | 0 / **3** | 2,171 / 2,399 | 0 / 0 |
| CROSS | 2,250 / 2,494 | 0 / **5** | 2,237 / 2,508 | 0 / **2** |
| OUTSIDE | 1,683 / 2,289 | 435 / 113 | 1,622 / 2,335 | 500 / 134 |
| REMOVED | 1,799 / 2,149 | 16 / 5 | 1,765 / 2,175 | 32 / 9 |
| BRANCH | 2,189 / 2,597 | 208 / 136 | 2,171 / 2,578 | 217 / 161 |
| SILENT | 1,415 / 1,695 | 0 / **2** | 1,339 / 1,695 | 0 / 0 |

- **B1 holds in all 12 worlds,** by +191 to +713.
- **B2 fails in four worlds,** where the current method makes no error: A119b INTACT, CROSS and SILENT, and public CROSS.
- **Round 10's rule and copy level give the same verdict, failing the same cells.**
- **No single family carries the failure.** Rescored with each family left out in turn, A119b fails every time; without
  DUX4 its errors are INTACT 1, CROSS 3, SILENT 1 against 0.

**The errors are O2 assignments at margin 2 to a near-identical wrong copy, never O3:**
- **DUX4 (80 copies):** a read's sequencing edits coincide with another unit's differences.
  - A DUX4.9 read is 2 edits from DUX4.19 and 4 from its own copy.
  - A DUX4.2 read matches DUX4.22 exactly and its own copy at 2.
  - With 80 near-identical competitors, chance coincidences beat a fixed margin of 2.
- **KRTAP2 (INTACT):** a KRTAP2.0 read aligns perfectly to the genome but extends past KRTAP2.0's annotated end, and only
  KRTAP2.1's longer isoform covers it. Its path distances are 2 against 85: structure decides, not sequence.
- **CDRT15 (CROSS):** CDRT15.0 and CDRT15.2 reads express structures annotated only on CDRT15.1 and CDRT15.3. The union paths
  do not give them those structures, as in round 10's WASH case.
- **The current method abstains on all of them.**

**Reported beside the bars:**
- **Arm I:**
  - **REMOVED:** all 4 detectable removed copies called under both profiles (DUX4.1, HLA-DRB.1, KRTAP1.1, KRTAP2.1), and none
    of the 4 within delta.
  - **BRANCH:** novel-exon copies of CDRT15, ETD, HLA-DRB, KRTAP2 and NAIP called under both profiles.
  - **OUTSIDE:** 3 (A119b) and 2 (public) PSV calls in KRTAP2 on present or uncatalogued reads (55 and 48 reads), where B2 still
    holds.
- **O1:**
  - **SILENT:** 6 of 33 silent copies dropped; 2-4 kept as undecidable (DUX4, ETD, NAIP, OR11H);
  - **OUTSIDE:** 5 new copies (DUX4, CDRT15, ETD, HLA-DRB, KRTAP2), two of them twins.
- **Cost: v5 is independent of the current method but not cheaper on an 80-copy array.**
  - **The integrated steps, summed over six worlds:** 1,095 s (A119b) and 1,624 s (public).
    - O3: 791 and 1,355 s, with one world's O3 alone over 590 s; it hosts reads on all 80 seeded units, and each host aligns
      the family pool.
    - O2: 174-208 s. O1: 90 s.
  - **The current method:** 383 and 356 s, peaking at 6.8-7.0 GB per world.
  - **Wall:** 1,098 s and 1,257 s per profile.

### Independent review (2026-10-16/17, from 23:32) and corrections

Fresh subagent; full report `/mnt/linuxdisk/tmp/o3_graph/review12/REVIEW.md`.

**The FAIL stands.**
- The reviewer's own scorer reproduces every bar, cell for cell, under all 3 rules.
- 00c241db reproduces O1, O2 and the fallback byte-identically.
- **The protocol is clean apart from Amendment 1, which could not change anything.** The crashed public run and its repeat
  are byte-identical in 71 of 72 outputs; the 72nd was never written by the crash.
- **Amendment 1's timing:** it was committed at 23:09:51, not at its 23:12 heading. A119b's bars and its error reads had been
  read first, because one tool output delivered both.
- **v5 uses nothing from the current method** (confirmed).

**Corrections:**
1. **The failure's mechanism was wrong for 5 of its 7 distinct reads (Critical).** The 12 failing error instances are 7 reads,
   repeated across worlds that share their reads.
   - **DUX4 (2 reads):** margin-2 coincidences, as stated. The genome's own primary alignment also puts both reads on the
     wrong unit (MAPQ 28-29).
   - **KRTAP2 (1 read):** it does not run past an annotated end.
     - It is KRTAP2.0's isoform with its intron RETAINED, DNA-identical at KRTAP2.0 and KRTAP2.1 (two NM-0 alignments), so
       abstaining is the only right answer.
     - Its 2 vs 85: KRTAP2.0's spliced path clips 83 bases, and a clip over 50 costs its length. KRTAP2.1's long-exon
       isoform absorbs the read's tandem-repeat difference as one 15-bp insertion, which costs one event.
   - **CDRT15, A119b (2 reads):** retained introns again, not missing structures.
     - The union does give CDRT15.0 the read's structure. The paralog's path absorbs the 265-bp retained intron as ONE
       142-bp insertion (gap-compressed 9), while CDRT15.0's own paths clip it at full cost.
     - One read went to a sole candidate: an O1-observed chain on CDRT15.1, learned from CROSS reads.
     - **The error floor contributes:** both reads sit at 0.0100 per base, between delta and delta + floor. Under delta
       alone they would be `unexplained` and placed correctly by the genome fallback.
   - **CDRT15, public (2 reads):** structures absent from the union, as stated.
   - **The current method abstains on the DUX4 and KRTAP2 reads, but places all 4 CDRT15 reads correctly.**
   - Retained introns cause 3 of A119b's 5 error reads. A119b retains introns at 6.4x the public rate.
2. **Leave one family out (Important):** public fails only through CDRT15. Without HLA-DRB, BRANCH B2 fails under both
   profiles.
3. **DUX4's novel-exon copy is missed, unreported (Important).** O2 places 57 and 62 of its reads that CARRY the novel exon,
   because a 150-bp novel exon costs one gap-compressed edit.
4. **The cost mechanism was wrong (Important).**
   - GraphAligner aligned 4% of the DUX4 reads in public OUTSIDE (61-96% elsewhere). 3.44 Mbp of whole reads became "novel
     segments".
   - The branch arm's all-against-all minimap2 then took 553 s and 1.7 GB, producing a 3,974-read call that the genome check
     removed.
   - O3's time tracks this novel-segment volume, not the 80 hosts, and peaked at 6.15 GB under the light lock.
5. **Arm I's KRTAP2.1 "detection" in REMOVED is identical to its OUTSIDE false calls (Important).**
   - KRTAP2.1's 52-57-bp exons in 15-bp tandem repeats mis-align even when the copy is present (INTACT median `de` 0.0137).
   - The OUTSIDE call members are all uncatalogued KRTAP2.1 reads: 35 and 24 distinct reads, not 55 and 48.
6. **The zero O3 errors in INTACT, CROSS and SILENT rest on the genome check alone (Important).** Before it, arm I's calls
   hold 1,203-1,462 present-copy reads per world, mostly DUX4 branch clusters of reads GraphAligner could not align.
7. **Minor:**
   - **O1 seeded a false new DUX4 copy in REMOVED** from the absent copy's reads (harmless to the bars).
   - **KRTAP1.1, NAIP.0 and OR11H.2 are unseeded outside copies.**
   - **The two twins keep few of their own reads.**
   - **The current method peaked at 6.2-6.5 GB under public**, not 6.8-7.0.
   - **The fates per read class are tabulated only in the review**, for the third round running.

**What this changes:**
- **Gap compression conflates sequencing indels with structure.** A retained intron or a novel exon becomes one edit when
  minimap2 bridges it, while the same sequence costs its full length when minimap2 clips it. That asymmetry, not near-identity,
  decided 5 of the 7 error reads, and it hides novel exons from O3.
- **The branch arm has no defence against GraphAligner failing on a graph;** its all-against-all then sets O3's cost.


---

## PREREG_o3_graph_round5_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round5_2026-10-16.md`)*

## Graph O2 with a structure-free fallback, and cross-family-guarded O3, on a sixth handpicked family set (2026-10-16)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal:** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to improve the graph
based O2 and O3". Round 4 (`docs/PREREG_o3_graph_round4_2026-10-15.md`, register 1302) claimed O3 for families without
sister-family hosting, and did not claim O2. Both failures had one root: the graph method gave up on reads it could not
explain within a family.
- **O2 failed A1 (0.969):** in CROSS it abstained on paralog-structure reads it does not thread, which the current O2 places
  correctly.
- **O3:** on development set 2 it pooled another family's reads (NPIPA.5 hosted on NPIPB) into a family's unexplained set
  and called them.

### Methods under test

**O2: graph O2 with a fallback** (`o2graph.py run WORLD unionf 2`, after `unionp 2` and `seqsupport`), `fallback_call`:
- padded union O2 decides first (round 4);
- a read it cannot explain in this family is `other_family` when another family's union explains it (within delta);
- otherwise the read is decided on its spliced alignments to each catalogued copy's genomic locus (`locus`: any canonical
  junction, the same margin and delta).

**O3: K** (`detect.py`): as round 4's H, but its pool is the reads the fallback O2 leaves `unexplained`. So another family's
reads and reads the family's own loci explain stay out.

Beside them: the current O2 and CUR; padded union O2 and H (round 4); locus O2; S23, S3, S123m; Kn (K before the genome
check).

### Worlds

As round 4, seed 61: INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon), and CROSS
(`--projection liftover`).

### Families (held out, picked before any run)

The candidates were taken in a fixed order: the TRIM and LILR subfamilies first, so that sister families could occur, then
SPDYE, STAG3L, MT1, GSTM, CEACAM, FAM230, HERC2P, CYP2D, BTN3A, SPANXN, GUSBP, ANKRD20A, FRG2, GOLGA6 and PMS2. Too few
qualified at first, so the last five were added in order, still before any run. The set is the first 8 with >= 3 copies.

| family | copies / isoforms |
|---|---|
| TRIM64 (64, 64B, 64C) | 3 / 3 |
| TRIM51 (TRIM51 and 8 pseudogenes, chr11 and chr2) | 9 / 12 |
| SPDYE | 10 / 17 |
| STAG3L | 3 / 12 |
| FAM230 | 9 / 20 |
| HERC2P | 3 / 4 |
| GUSBP | 3 / 3 |
| ANKRD20A (8 chromosomes) | 10 / 12 |

`/mnt/linuxdisk/tmp/o3_graph/families5.json` (the earlier sets identical to `families4.json`).

### Known before running (development: the five earlier family sets)

O2, the fallback graph O2 against round 4's bars (re-labelled C1-C4 below):

| set | C1: lowest world ratio of correct calls to current O2's | C2: lowest precision | C3: REMOVED and BRANCH fewer; pooled non-catalogued ratio | C4: vs locus O2 | correct: fallback / padded union / current |
|---|---|---|---|---|---|
| dev | 1.025 | 1.0 | yes; 0.52 | x1.049 | 11,389 / 11,315 / 10,939 |
| dev-2 | 1.083 | 1.0 | yes; 0.41 | x1.126 | 12,702 / 12,701 / 11,530 |
| dev-3 | 1.065 | 1.0 | yes; 0.45 | x1.055 | 7,115 / 7,106 / 6,678 |
| dev-4 | 0.997 | 1.0 | yes; 0.54 | x1.056 | 9,220 / 9,185 / 9,230 |
| dev-5 (round 4's set) | 1.002 | 1.0 | yes; 0.29 | x1.059 | 8,928 / 8,814 / 8,906 |

- Every paralog-structure read the fallback assigns is right: 957, 612, 270, 222 and 471 of them; there are no wrong-family
  calls.
- The fallback places more non-catalogued reads than padded union (dev-4: 479 vs 375). C3's pooled ratio rises towards 0.54.

O3, K on dev-2, the only development set with sister-family hosting (550-790 reads per world hosted on another family):

| world | H false-call families | K false-call families | K detected |
|---|---|---|---|
| INTACT | 1 (NPIPB) | 0 | - |
| CROSS | 1 | 0 | - |
| OUTSIDE | 1 | 0 | - |
| REMOVED | 1 | 0 | 2 / 2 detectable |
| BRANCH | 0 | 0 | 8 / 8 |

K on the other development sets was still running at commit; it is reported under Results.

### Bars (held out: the 8 families above, seed 61; wrong-family calls count as wrong)

| id | bar |
|---|---|
| C1 | fallback graph O2 correct >= 0.98 x the current O2's, in every world |
| C2 | fallback graph O2 precision >= 0.995 over reads of catalogued copies, in every world |
| C3 | in REMOVED and in BRANCH it places fewer reads of non-catalogued copies than the current O2, and <= 0.6 x as many over OUTSIDE, REMOVED and BRANCH together |
| C4 | fallback graph O2 correct >= 1.03 x locus O2's, pooled over the five worlds |
| D1 | K detects >= 6 of the 8 BRANCH absent copies |
| D2 | K has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| D3 | K detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| D4 | K detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |
| D5 | in every family-world with >= 100 reads hosted from another family, K has no false call; not measurable without one |

Detectable (D3): the removed copy's reads have a median primary `de` on the masked genome > delta (`summary.py o3`).

### Readings (fixed now)

1. **C1-C4 all pass:** claimed, "on simulated reads of handpicked real families, graph O2 that threads every isoform structure
   through every copy, and falls back to each copy's genomic locus for reads no structure explains:
   - keeps >= 98% of the current O2's correct assignments in every world, including when isoforms are annotated on another
     paralog;
   - has >= 99.5% precision;
   - places <= 60% as many reads of non-catalogued copies;
   - keeps >= 3% more correct assignments than a locus-only O2."

   Any failing bar: not claimed.
2. **D1-D4 all pass:** claimed, round 4's O3 sentence for K. D3 not measurable: claimed without the removed-copy clause. D3
   measurable and failing, or any other bar failing: not claimed. D5 is reported. If it is measurable and fails, the claim
   keeps round 4's scope (families without sister-family hosting).
3. Each claim stands or falls as a whole. A failed bar is reported by name.

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (held out, 2026-10-16)

Run: code and prereg `1c49489a`, the frozen copy_assign binary, seed 61. Outputs in `/mnt/linuxdisk/tmp/o3_graph/real/heldout5/`
(`r5_*.json`, `o3_score.json`). Summaries from `bench/o3_graph/summary.py` (`o2bars-unionf`, `o3`).

**O2** (correct / assigned; wrong-family calls counted):

| world | current | **fallback graph O2** | locus O2 | non-catalogued placed: current / fallback |
|---|---|---|---|---|
| INTACT | 2,357 | 2,546 / 2,546 (x1.080) | 2,456 | - |
| CROSS (liftover) | 2,793 | 3,005 / 3,006 (x1.076) | 2,922 | - |
| OUTSIDE | 2,066 | 2,260 / 2,260 (x1.094) | 2,237 | 199 / 449 |
| REMOVED | 2,066 | 2,260 / 2,260 (x1.094) | 2,237 | 233 / 449 |
| BRANCH | 2,357 | 2,546 / 2,546 (x1.080) | 2,456 | 104 / 0 |

- On this set the uncatalogued (OUTSIDE) and removed (REMOVED) special copies are all within delta of a catalogued paralog: the
  removed copies' median `de` is 0.0013-0.0077, so none is detectable.
- Their reads are explained by that paralog. The fallback, like any distance-based O2, places 449 of them on it (wrong by
  construction). The current O2 abstains on more of them: 199 and 233.

**O3** (8 families per world):

| | CUR | S23 | S3 | H | **K** | Kn (before the genome check) |
|---|---|---|---|---|---|---|
| BRANCH detected | 0 | 8 | 0 | 8 | **8** | 8 |
| INTACT + OUTSIDE + CROSS false-call families (of 24) | 2 | 0 | 0 | 0 | **0** | 5 |
| REMOVED + BRANCH false-call families | 2 | 0 | 0 | 0 | **0** | 4 |
| REMOVED detected (0 detectable) | - | - | - | - | - | - |

**Sister-family hosting (D5).** TRIM64 hosts 143-230 reads from another family (TRIM51) in every world, and TRIM51 hosts 60 in
REMOVED. K has no false call in any of them; Kn has TRIM64 or TRIM51 false calls in all five worlds, which the guard and the
genome check remove. H, without the guard, has none either, so on this set the genome check alone sufficed. On development
set 2 it did not: H had NPIPB false calls, and K none.

| bar | value | |
|---|---|---|
| C1 | lowest world ratio 1.076 | pass |
| C2 | lowest precision 0.9997 (CROSS 3,005 / 3,006) | pass |
| C3 | REMOVED 449 vs 233 (more); pooled 898 / 536 = 1.68 | **fail** |
| C4 | 12,617 vs 12,308: x1.025 | **fail** |
| D1 | 8 / 8 | pass |
| D2 | 0 of 24 | pass |
| D3 | 0 detectable removed copies | not measurable |
| D4 | K 8 = S23 8, false-call families 0 = 0 | pass |
| D5 | measurable (TRIM64, TRIM51); K no false call | pass |

#### Verdict (2026-10-16)

**O3: claimed without the removed-copy clause (reading 2: D1, D2, D4 pass; D3 not measurable).** "Graph O3 that clusters what
no family's graph, nor the family's own loci, explains detects novel-exon absent copies (8 of 8; the current method 0). It
makes no false call over intact, uncatalogued-paralog and paralog-annotation families (0 of 24; the current method 2),
including families that host another family's reads (D5)."

Together with round 4, which claimed the removed-copy clause for H (6 of 6) on families without sister-family hosting:
- the graph-based O3 holds on two fresh sets;
- K detects everything H detects on all five development sets, and has no false call on any of them.

Attribution, as in round 4: the novel exons come from the branch calls, and the false-call drop from the genome check plus,
for sister families, the cross-family guard.

**O2: not claimed (C3, C4 fail).** The fallback graph O2 assigns 7.6-9.4% more reads correctly than the current O2 in every
world, at precision >= 0.9997, including the paralog-annotated isoforms (C1, C2). But:
- where an uncatalogued copy is within delta of a catalogued one, it places that copy's reads on the paralog, as any
  distance-based O2 would, while the current O2 abstains on more of them (C3);
- it keeps only 2.5% more correct calls than locus O2 here (C4).

Across six fresh family sets graph O2 has traded one failure for another:
- round 4: abstaining on structures it does not thread;
- round 5: assigning reads of near-identical uncatalogued copies.

No version has met all of its pre-registered bars.

### Review corrections (round-5 review, 2026-10-16; post hoc)

A fresh reviewer re-derived every number above with independent code, and all of them reproduce. Timing is clean: commit
00:52:34.052, first held-out file 7 ms later. The family selection reproduces exactly (first 8 with >= 3 copies, in GROUPS
order, from `families5_all.json`). No code bug changes a number. The corrections are about what the held-out set showed.

**D5 passed without exercising the guard (I1).**
- None of the TRIM51 reads hosted on TRIM64 (148-230 per world) aligns to TRIM64's loci, so no sister-family PSV cluster forms
  there. Development set 2 had 30-43 per world.
- Hn's and Kn's sister-family calls are identical. All of them are branch calls, and the genome check drops all of them
  (median `de` 0.0011-0.0041, no novel bases).
- So "which the guard and the genome check remove" was wrong: the guard removed nothing.
- "Measurable (TRIM64, TRIM51)" was wrong too: TRIM51 hosts 60 foreign reads, under the 100 needed; only TRIM64 qualifies.
- D5 is passed but says nothing about the guard. The sister-family evidence for K is development set 2 only.

**The guard covers PSV clustering only (I2).** K's branch calls are built from all of a family's hosted reads, as S2's are, not
from K's pool. So another family's reads with >= 50 genome-unaligned bases, or `de` > delta, could still be called; only the
genome check stands in the way. The prereg's "another family's reads ... stay out" holds for K's PSV clustering only. Not
changed after the fact: it would be a method change for a new prereg.

**Wording, attribution and history (I3, I4).**
- The claim sentence should be the pre-registered one (below).
- On this set K, H and S23 are equal in every cell: 8 / 8 BRANCH, 0 / 24 false-call families, 0 in REMOVED + BRANCH. S3,
  the current method plus the genome check, is also 0 / 24.
- All 8 detections are branch calls that S23 also makes, and the removed-copy clause had no data. So none of K's own parts
  (the fallback pool, the guard, the mosaic rule, pooling) was tested on held-out data.
- On development data K adds over S23 only in REMOVED (GAGE, GOLGA8, PSG, CGB, KIR, PRAMEF), exactly the clause this set
  could not measure.
- This is round 5, and the detector changed after each round. The guard fixes false calls that round 4's mosaic extension
  introduced; S23 has none on dev-2.
- "Holds on two fresh sets" mixed H (round 4) and K (round 5).

**Where C3's failure comes from (I5).** Placements of the 481 non-catalogued reads per world:

| world | current O2 | padded union | locus O2 | fallback |
|---|---|---|---|---|
| OUTSIDE | 199 | 303 | 437 | 449 |
| REMOVED | 233 | 303 | 437 | 449 |

- The fallback adds 146 per world: its own `locus` calls on reads padded union left `unexplained` (TRIM51.7 60, HERC2P 59,
  STAG3L 27). It cannot tell a catalogued copy's unthreaded structure from a near-identical uncatalogued copy.
- Its own calls on held-out: 48 correct against 292 non-catalogued reads placed. On the development sets: 233 against 259.
- Counting non-catalogued placements, precision is 0.834 in OUTSIDE and REMOVED, against the current O2's 0.912 and 0.899.
  The ">= 0.9997" was C2's catalogued-only figure.
- On this set padded union alone would also fail C3 (REMOVED 303 > 233; pooled 1.13) and C4 (x1.021). It would pass C1
  (lowest 1.062, CROSS).
- So C3 fails mostly because of this family set, and the fallback makes it worse.

**The O2 history, by round (I6).** Five fresh held-out sets, each failing on a different bar:
- **Round 1:** path O2 passed its own bars (O2a-O2d). The review then found the gain annotation-decided, and Amendment 1's
  CROSS failed (precision 0.42 on paralog-structure reads).
- **Round 2:** union O2, U4 0.981.
- **Round 3:** padded union, no gain (x0.995).
- **Round 4:** padded union, A1 0.969 in CROSS.
- **Round 5:** the fallback, C3 and C4.

"Traded one failure for another" shows only on development data: the fallback raises CROSS (dev-5: 2,330 vs 2,252) and
raises non-catalogued placements (dev-4: 479 vs 375).

**Smaller.**
- The "dev" row pools six worlds (CROSS and CROSSL).
- Not fixed, recorded:
  - `summary.py` labels the fallback's output "padded union" and its bars A1-A4;
  - the guard does not stop a `locus` call when the read's own family's union also fails (none occurred);
  - there is no table-level test of `run_fallback` (the reviewer re-derived every row);
  - the selection step is not a committed command (it reproduces exactly from the rule above);
  - a future prereg should require >= 3 detectable removed copies before running, since detectability is computable
    without reads.

#### Verdict, revised after the review (2026-10-16)

**O3: claimed without the removed-copy clause (reading 2), in the pre-registered words.** "On simulated reads of handpicked real
families, graph O3 that clusters a family's O2-unexplained reads together, with the mosaic rule, detects novel-exon absent
copies with at most one false call over intact, uncatalogued-paralog and paralog-annotation families." Observed: 8 of 8 and
0 of 24.

What this set tested: the branch calls (novel exons) and the genome check (no false call), exactly as S23 and S3 do. It did not
test K's own parts. The sister-family evidence for the guard is development set 2 only (K 0 false calls, H 4).

**O2: not claimed (C3, C4 fail).**
- The fallback assigns 7.6-9.4% more reads correctly per world, at precision >= 0.9997 over reads of catalogued copies.
- But its own `locus` calls are mostly reads of near-identical uncatalogued copies. Counting those, its precision is 0.834,
  against the current O2's 0.90-0.91 in OUTSIDE and REMOVED.
- Across five fresh held-out sets no graph O2 version has met all of its pre-registered bars.

### Post hoc: can any operating point of graph O2 dominate the current O2? (2026-10-16; all sets seen, no claim)

Worlds with an uncatalogued copy (OUTSIDE, REMOVED); a placement of a non-catalogued copy's read counts as wrong. Correct /
wrong:

| set | world | current O2 | fallback graph O2 |
|---|---|---|---|
| round 2's | REMOVED / OUTSIDE | 1,034 / 162; 1,034 / 119 | 1,102 / 119 (dominates both) |
| round 3's | REMOVED / OUTSIDE | 1,624 / 373; 1,621 / 241 | 1,619 / 235 (5 fewer correct, far fewer wrong) |
| round 4's | REMOVED / OUTSIDE | 1,437 / 121; 1,434 / 60 | 1,440 / 60 (dominates both) |
| round 5's | REMOVED / OUTSIDE | 2,066 / 233; 2,066 / 199 | 2,260 / 449 (trade-off) |

Could a per-read error test separate a near-identical uncatalogued copy's reads from sequencing error? Reads assigned to their
own copy sit at a median 1.6-1.7 differences per kb on the copy's locus, with a heavy tail (95th percentile 6-14). The
uncatalogued copies' reads sit at 4-5 per kb. No cutoff from 0.5 to 9.5 per kb makes graph O2 dominate where it does not
already.

So graph O2 is at least as good as the current O2 on most sets and worlds, and better in precision. It loses precision exactly
where an uncatalogued copy is within sequencing-error range of a catalogued one (round 5's set), and nothing tried here tells
those reads apart.


---

## PREREG_o3_graph_round6_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round6_2026-10-16.md`)*

## O2 + an unannotated-copy flag, with heterozygous alleles, on a seventh handpicked family set (2026-10-16)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User question (2026-10-16):** "can we check if O2 is doing that because it is detecting unannotated copies? In that case that
could be useful". The answer was to pre-register the flag on a fresh family set, with alleles so that the copy-versus-allele
cost is measured; the user agreed ("yes please").

### Background (post hoc, round 5's set and earlier, no alleles)

O2's wrong placements in OUTSIDE and REMOVED worlds are reads of the uncatalogued copy, placed on its nearest catalogued
paralog. As a group they often form a pure PSV cluster whose copy-private alleles are linked. That is an unannotated copy's
fingerprint, usually at 0.1-0.5% divergence, below delta (0.96%).

Across six family sets, the flag defined below:
- covered 24 of 48 removed copies, 13 of them below delta;
- held 45% of O2's misplaced reads;
- falsely flagged 4 of 48 intact families.

There are no alleles in those worlds. In real data a heterozygous site gives the same signature (linked differences shared by
about half of a copy's reads). This round measures that cost.

### Simulation: heterozygous alleles (`realworld.py --het RATE`)

The donor is diploid. Each copy (catalogued or not) has heterozygous sites at RATE of its exonic bases, and its haplotype 2
carries the alternative bases. Each read comes from haplotype 1 or 2 with probability 1/2 (truth column `hap`). **RATE = 0.001**,
human heterozygosity. With `--het 0` a world is byte-identical to the earlier simulator's.

### Methods under test

- **The flag** (`o2graph.flag_reads`): a PSV cluster of the current method's per-host clustering (`copy_assign
  --psv-clusters`) whose copy-private alleles are linked (`unexplained`), with >= 20 reads, whatever its divergence.
- **O2 + flag** (`o2graph.py run WORLD flag unionp`): padded union O2 (round 3), with every read of a flagged cluster reported
  as an unannotated-copy candidate instead of being assigned.

Padded union is the base because, without alleles, padded union + flag had at least the current O2's correct calls (within
0.4%) and fewer wrong ones in all 18 development worlds.

### Worlds

INTACT, REMOVED (special copy hard-masked, own index) and OUTSIDE, all with `--het 0.001`, seed 71.

### Families (held out, picked before any run)

The candidates were taken in a fixed order: FRG2, GOLGA6, PMS2, CCL3/4, DEFA, BCRP, GGT, CFHR, PRB, ANKRD30, FRG1, GSTT2, CTAGE,
RIMBP3. FRG2, GOLGA6 and PMS2 were round-5 candidates past its first 8 and never built into a world. The set is the first 8
with >= 3 copies:

| family | copies / isoforms |
|---|---|
| FRG2 | 5 / 9 |
| GOLGA6 (6B, 6C, 6D) | 3 / 5 |
| PMS2 (PMS2 + 9 pseudogenes) | 10 / 20 |
| CCL34 (CCL3, CCL3L1, CCL3L3) | 3 / 9 |
| DEFA (DEFA1, DEFA1B, DEFA3) | 3 / 4 |
| BCRP | 5 / 5 |
| GGT | 5 / 10 |
| FRG1 | 10 / 10 |

`/mnt/linuxdisk/tmp/o3_graph/families6.json` (the earlier sets identical to `families5.json`).

### Known before running (development: round 5's families with `--het 0.001`, seed 61; `real/dev_het/`)

| | INTACT | REMOVED | OUTSIDE |
|---|---|---|---|
| flagged clusters: copy / allele / other | 0 / 13 / 8 | 6 / 10 / 7 | 6 / 10 / 7 |
| uncatalogued copies flagged | - | 2 of 8 | 2 of 8 |
| families with a false flag | 7 of 8 | 7 of 8 | 7 of 8 |
| correct calls: current O2 / padded union / padded union + flag | 2,356 / 2,542 / 1,886 | 2,065 / 2,256 / 1,727 | 2,065 / 2,256 / 1,727 |
| wrong calls (incl. non-catalogued): current / padded union / + flag | 9 / 0 / 0 | 226 / 297 / 165 | 192 / 297 / 165 |

With alleles the flag fires on one haplotype of most catalogued copies. It costs 26% of padded union's correct calls in
INTACT, and the uncatalogued copies' reads split by haplotype, so fewer copies are covered. The same families without alleles
had 2 false-flag families and 4 of 8 copies flagged.

### Bars (held out: the 8 families above, seed 71, `--het 0.001`)

| id | bar |
|---|---|
| H1 | the flag covers >= 50% of the uncatalogued copies (REMOVED and OUTSIDE together, 16): >= 80% of a copy's hosted reads in flagged clusters that are >= 90% its reads |
| H2 | the flag falsely flags (a cluster not >= 90% one uncatalogued copy's reads) <= 1 of the 8 INTACT families |
| H3 | O2 + flag keeps >= 0.98 x padded union O2's correct calls in INTACT |
| H4 | O2 + flag places fewer wrong reads (non-catalogued + wrong copy + wrong family) than the current O2 in REMOVED and in OUTSIDE |

Truth classes come from `summary.py flags`: `copy`, `allele` (>= 90% one catalogued copy's reads, >= 90% of those from one
haplotype) and `other`.

### Readings (fixed now)

1. **H1-H4 all pass:** claimed, "with realistic heterozygosity, flagging PSV clusters of linked copy-private alleles below
   delta turns O2's misplaced reads of unannotated copies into candidates. It covers >= half of them, falsely flags <= 1 of 8
   intact families, and keeps >= 98% of O2's correct calls."
2. **H2 or H3 fails:** not claimed. "Without allele information (DNA of the same individual, or phasing), the flag cannot
   tell an unannotated copy below delta from a heterozygous allele." The development data predicts this reading.
3. **H1 fails alone:** not claimed; the flag misses unannotated copies even where alleles are not the problem.

Scope: simulated reads (curated CAT isoforms, HiFi errors, heterozygous sites at 0.1%) on the real CHM13 genome. No allelic
expression imbalance, no somatic or editing variation.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (held out, 2026-10-16)

Run: code and prereg `ef78d6ed`, the frozen copy_assign binary, seed 71, `--het 0.001` (60 heterozygous sites). Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout6/` (`flags.json`, `o2_*.json`); summaries from `summary.py flags` and `o2graph.py
score`.

| | INTACT | REMOVED | OUTSIDE |
|---|---|---|---|
| flagged clusters: copy / allele / other | 0 / 11 / 2 | 5 / 7 / 2 | 5 / 7 / 2 |
| uncatalogued copies flagged | - | 2 of 8 | 2 of 8 |
| families with a false flag | 4 of 8 (FRG2, GGT, GOLGA6, PMS2) | 4 | 4 |
| correct / wrong: current O2 | 1,983 / 0 | 1,683 / 157 | 1,683 / 157 |
| correct / wrong: padded union O2 | 2,104 / 0 | 1,790 / 118 | 1,790 / 118 |
| correct / wrong: padded union + flag | 1,701 / 0 | 1,505 / 117 | 1,505 / 117 |

"Wrong" counts reads of non-catalogued copies placed, wrong-copy and wrong-family calls.

| bar | value | |
|---|---|---|
| H1 | 4 of 16 uncatalogued copies covered (0.25) | **fail** |
| H2 | 4 of 8 INTACT families falsely flagged | **fail** |
| H3 | 1,701 / 2,104 = 0.81 of padded union's correct calls kept | **fail** |
| H4 | 117 < 157 in REMOVED and in OUTSIDE | pass |

#### Verdict (2026-10-16)

**Not claimed (reading 2: H2 and H3 fail).** Without allele information (DNA of the same individual, or phasing), the flag
cannot tell an unannotated copy below delta from a heterozygous allele.
- At human heterozygosity (0.1%), 11 of the 13 clusters it flags in intact families are one haplotype of a catalogued copy.
- It costs 19% of O2's correct calls.
- The unannotated copies' own reads split by haplotype, so it covers only 4 of 16.

So the answer to the question that started this round: O2's misplaced reads do carry an unannotated-copy signature, but at
the divergence where O2 misplaces them, alleles carry the same one.

**Reported beside the bars (no claim).** On this fresh set, padded union O2 without the flag beats the current O2 on both
counts in every world:
- INTACT 2,104 vs 1,983 correct (no wrong);
- REMOVED and OUTSIDE 1,790 vs 1,683 correct, and 118 vs 157 wrong.

Elsewhere, from the round-6 development tables (corrected after a first draft overstated it):
- it dominates in every world on the original development set and on round 2's set;
- on round 1's set it dominates except in OUTSIDE (87 vs 81 wrong);
- on round 3's and round 4's sets it is within 0.4% of the current O2's correct calls with fewer or equal wrong ones;
- on round 5's set, with or without alleles, it places more wrong reads (its uncatalogued copies are near-identical to
  catalogued ones). That is a candidate O2 claim for its own
pre-registration on a new family set: "padded union O2 dominates the current O2". It was not this round's test, so it is not
claimed here.

### Review corrections (round-6 review, 2026-10-16; post hoc)

A fresh reviewer re-derived every number above with independent code; all reproduce. Timing is clean (commit 09:32:19.491,
`heldout6/` 18 ms later; nothing above "## Results" changed). No Critical finding. The corrections are about why the bars
failed and what the flag did.

**Why H1 fails (I1).** It was not mainly haplotype splitting:
- Five of the six copies not covered have no linked copy-private alleles at all, so with MIN_ALT = 2 they can never be
  flagged, with or without alleles:
  - BCRP6: 6 alts, none unexplained;
  - CCL3L1: one unexplained alt;
  - DEFA1B: one private site, and DEFA got no het site;
  - FRG1.2 and PMS2.8: their haplotype-1 clusters are explained.
- Only GOLGA6B is a near miss: 47 / 60 = 78%, split over two hosts.
- GGT2 did split by haplotype (32 + 27) and was still covered.

So H1 also fails for reading 3's reason (copies with no PSV signature).

**The flag contributed almost nothing to O2 here (I2).**
- Padded union misplaces 118 reads: FRG1.2 60, CCL3L1 57, GOLGA6B 1. Only 1 of them is in a flagged cluster, against 45% on
  the six earlier sets.
- Of the 286 reads labelled `flagged` in REMOVED and OUTSIDE, 285 were correct calls of catalogued copies.
- So H4 passes because of padded union, not the flag.

**Dominance scope (I3).**
- The side observation holds for INTACT, REMOVED and OUTSIDE only. In CROSS, padded union keeps 0.983 (round 3's set) and
  0.969 (round 4's: that round's A1 failure) of the current O2's correct calls; round 4's REMOVED is 0.9958.
- This set has no CROSS world, and its REMOVED and OUTSIDE are identical for O2 (2,644 reads, same hosting).
- A dominance prereg needs a CROSS world.

**Realism (I4).** The verdict's direction is robust, but "19%", "4 of 8" and "11 of 13" depend on these design choices:
- Haplotype 1 is the reference exactly, and every alt allele sits on haplotype 2. A real donor's homozygous differences from
  the reference would make whole copies flaggable, so the true cost is probably higher. Random phase would lower the flag
  count; the reviewer's estimate is H2 still about 4 of 8 and H3 about 0.87, both failing.
- Het sites are a deterministic count (round(0.001 x exonic length)). Copies under 1.5 kb of exon get at most one site and can
  never be flagged, and the CCL3 copies sit just under. A Poisson draw would give about 6 of 8 families with a flaggable copy.
- Expression is uniform: about 30 reads per haplotype, above the 20-read minimum.

**Smaller.**
- **M1.** With het 0 every product is byte-identical except the truth file, which gains a `hap` column (all 1).
- **M2.** `classify_cluster`'s allele rule counts reads that cover no het site, and labels a >= 90% haplotype-1 cluster as
  an allele. By the truth, all 13 INTACT flags here are allele clusters (two were labelled `other`).
- **M3.** "13 below delta" in the background counts by PSV est_div; by O3's genome-`de` measure it is 10.
- **M4.** Pre-run statements:
  - "fewer wrong in all 18 development worlds": 5 tie;
  - "within 0.4%": the worst is 0.42%;
  - "fires on most catalogued copies": 13 of 44 held-out copies, which is most of those with >= 2 het sites.
- **M5.** `apply_flag` labels only reads O2 had assigned. Members of flagged `copy` clusters that O2 abstained on are not
  listed as candidates.
- **M6.** REMOVED and OUTSIDE are identical here, so "4 of 16" is 2 of 8 distinct copies counted twice.

#### Verdict, revised after the review (2026-10-16)

**Not claimed (reading 2: H2 and H3 fail; H1 also fails for reading 3's reason).**

The answer to the question that started this round:
- On the six earlier sets, about 45% of O2's misplaced reads carried an unannotated-copy signature.
- On this fresh set, O2's misplacements come from copies with no PSV signature at all, while the flag fires on heterozygous
  alleles: all 13 intact flags are allele clusters.
- So the flag is not usable without allele information (DNA of the same individual, or phasing), and it is not what explains
  O2's misplacements.


---

## PREREG_o3_graph_round7_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round7_2026-10-16.md`)*

## Integrated O1 -> O2 -> O3 on an eighth handpicked family set, with alleles (2026-10-16)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User question (2026-10-16):** "would fully integrating O1, O2 and O3 help? I think O1 might reduce the possible universe by
only seeding what is supported by RNA reads". The plan was approved: "yes, please do all of that".

**Why it might help (from rounds 1-6).** O2's errors come from three places:
- a copy present in the genome but missing from the catalog (OUTSIDE). Seeding from RNA support could add it: 48 of 56 such
  copies over seven earlier sets have most of their reads' primary alignments on their own locus;
- isoforms annotated only on a paralog (CROSS). A copy's observed structures could add them;
- copies absent from the reference, near-identical to a paralog, and alleles. Integration cannot fix these.

Shrinking the universe to RNA-supported copies could reduce ties. Every catalogued copy was expressed in the earlier worlds,
so this is measured here for the first time (SILENT).

### Methods under test

**O1 seeding** (`o1seed.py WORLD` -> `graph/seeded.json`), per family:
- **family reads:** reads with any alignment (primary or secondary) on a curated copy's locus (+-1 kb);
- **kept:** a curated copy is kept when >= 5 family reads have their primary alignment on it; otherwise it is DROPPED
  (silent). Every primary counts here, so an identical expressed copy, whose reads split at mapq 0, is never dropped;
- **new:** a locus off every family's catalog becomes a new copy when >= 5 family reads have their primary alignment there
  UNIQUELY (mapq >= 1). Adding a copy needs positive evidence; a mapq-0 primary is a coin flip;
- **structures:** each copy's curated isoforms, plus the intron chains seen in >= 3 near-uniquely placed (mapq >= 50) primary
  alignments that this copy overlaps most and that lie within its span (+-100 bp). A new copy has its observed chains only,
  from uniquely placed reads.

**Integrated O2** (`o2graph.py run WORLD integ 2`): padded union O2 (round 3: every structure threaded through every copy, path
ends padded with each copy's genome) on the seeded catalog, over the family reads.

**Integrated O3** (`detect.py`, arm **I**): as round 4's H (pooled, mosaic rule, branch calls, genome check), on the reads the
integrated O2 leaves `unexplained`.

Beside them: the current O2 and CUR; padded union O2 on the curated catalog; H and S23.

### Worlds

All with heterozygous alleles (`--het 0.001`), seed 81:
- INTACT;
- CROSS (`--projection liftover`);
- OUTSIDE;
- REMOVED (special copy hard-masked, own index);
- **SILENT**: a seeded third of each family's copies express nothing; all of them stay in the catalog.

**Seed (detectability pre-check, round-5 review M6).** `precheck.py` aligns each candidate special copy's main isoform to the
genome; its divergence is that of its best alignment away from its own locus. Seed 81 is the smallest seed >= 81 with >= 3
detectable families: SPRR2 0.067, TUBB8 0.086, WASH 0.0116. The others are within delta.

### Families (held out, picked before any run)

The candidates were taken in a fixed order: CTAGE, AQP7, PRSS, SPRR2, FAM74A, TUBB8, DDX11L, WASH, SIGLEC, IFITM, FAM138,
FAM157. CTAGE was a round-6 candidate past its first 8 and was never built into a world. The set is the first 8 with >= 3
copies:

| family | copies / isoforms |
|---|---|
| CTAGE | 4 / 4 |
| AQP7 (AQP7 + pseudogenes) | 4 / 7 |
| SPRR2 | 5 / 9 |
| FAM74A | 4 / 4 |
| TUBB8 (TUBB8 + pseudogenes) | 6 / 11 |
| DDX11L | 10 / 13 |
| WASH | 9 / 22 |
| FAM138 | 8 / 17 |

`/mnt/linuxdisk/tmp/o3_graph/families7.json` (the earlier sets identical to `families6.json`).

### Scoring

- A call is correct when it names the read's own copy. A new seeded copy is mapped to the truth copy whose locus it overlaps
  (scoring only).
- **Reads of an absent copy that the genome explains within delta** (undetectable by O3's definition: its reads' median `de`
  on the masked genome <= delta) cannot be told from the locus that explains them. Placing them there is counted apart
  (`unknowable_placed`), for every method alike, and reported; it is not an error. All other wrong calls are errors.

### Known before running (development)

Development: round 6's families with alleles (INTACT, REMOVED, OUTSIDE; plus new SILENT and CROSS worlds, seed 71), round 5's
families with alleles (`dev_het`), and round 4's families without alleles. Correct / wrong calls (unknowable placements in
brackets):

| set | world | current O2 | padded union O2 | **integrated O2** | ratio of correct calls to the current O2's | uncatalogued reads on their own copy | new copies (on a truth copy) |
|---|---|---|---|---|---|---|---|
| round 6 | SILENT | 1,380 / 0 | 1,509 / 0 | **1,569 / 0** | x1.137 | - | 0 |
| round 6 | CROSS | 2,131 / 0 | 2,318 / 0 | **2,315 / 0** | x1.086 | - | 0 |
| round 6 | INTACT | 1,983 / 0 | 2,104 / 0 | **2,104 / 0** | x1.061 | - | 0 |
| round 6 | REMOVED | 1,683 / 0 (157) | 1,790 / 0 (118) | **1,790 / 0 (165)** | x1.064 | - | 2 (0) |
| round 6 | OUTSIDE | 1,683 / 157 | 1,790 / 118 | **2,074 / 0** | x1.232 | 335 / 480 | 8 (8) |
| round 5 | INTACT | 2,356 / 9 | 2,542 / 0 | **2,542 / 0** | x1.079 | - | 0 |
| round 5 | REMOVED | 2,065 / 0 (226) | 2,256 / 0 (297) | **2,265 / 0 (323)** | x1.097 | - | 1 (0) |
| round 5 | OUTSIDE | 2,065 / 192 | 2,256 / 297 | **2,542 / 0** | x1.231 | 347 / 481 | 7 (7) |
| round 4 | INTACT | 1,855 / 0 | 1,850 / 0 | **1,850 / 0** | x0.997 | - | 0 |
| round 4 | CROSS | 2,325 / 0 | 2,252 / 0 | **2,313 / 0** | x0.995 | - | 0 |
| round 4 | OUTSIDE | 1,434 / 60 | 1,431 / 60 | **1,795 / 0** | x1.252 | 364 / 481 | 7 (7) |
| round 4 | REMOVED | 1,437 / 61 (60) | 1,431 / 0 (60) | **1,431 / 0 (60)** | x0.996 | - | 1 (0) |

O1 dropped all 12 silent copies and none of the 344 expressed copies over these worlds.

Found and fixed on development, before this commit:
- **New-locus support counts unique placements only (mapq > 0).** The first version seeded 10 genomic paralogs in round 6's
  CROSS from mapq-0 coin flips. They competed with the true copies: 2,253 correct against padded union's 2,318.
- **Each alignment counts for the one curated copy it overlaps most, and a structure only from alignments inside that copy's
  span.** CSAG's copies 1 and 2 are 91 bp apart, so their padded windows took each other's reads and structures, and 105
  reads tied.
- **A curated copy's observed structures need near-unique reads (mapq >= 50),** the user's suggestion against fragments of
  ambiguous placement.
- **A new locus still counts any unique placement (mapq >= 1).** With >= 50, real uncatalogued copies near-identical to a
  paralog (ANKRD20A.4, SPDYE.2, TRIM64.2, DEFA.1, FRG1.2, CSAG.3: all reads unique at mapq 1-49) were not seeded, and
  round 5's OUTSIDE recovered only 161 of 481 reads.
- **New copies in REMOVED (1-2 per world) are genomic paralogs that absorbed a removed copy's reads.** All of those copies
  are within delta of the genome, so the placements are counted as unknowable.

O3 (arm I) on these worlds was re-run with this final seeding after the commit; it is reported under Results. With the
previous seeding (unique placements for new loci, no ownership rule), arm I had:
- no false call in any of the 12 worlds;
- detections of 3 of 3 detectable removed copies in round 6's REMOVED and 6 of 6 in round 4's.

### Bars (held out: the 8 families above, seed 81, five worlds)

| id | bar |
|---|---|
| J1 | the integrated O2 has >= 0.98 x the current O2's correct calls, in every world (as rounds 4-5's A1 / C1) |
| J2 | the integrated O2 has <= the current O2's wrong calls, in every world |
| J3 | in OUTSIDE the integrated O2 assigns >= 50% of the uncatalogued in-genome copies' reads to their own (seeded) copy |
| J4 | O1 drops >= 80% of SILENT's silent copies and <= 10% of expressed copies (pooled over the five worlds) |
| J5 | in SILENT the integrated O2 has >= padded union O2's correct calls on the curated catalog |
| J6 | arm I has a false call in <= 1 of the 32 INTACT, CROSS, OUTSIDE and SILENT families |
| J7 | arm I detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |

### Readings (fixed now)

1. **J1-J5 all pass:** claimed, "on simulated reads of handpicked real families with heterozygous alleles, an O1 that seeds
   the copy catalog from RNA support makes O2 as right (>= 98% of the correct calls) and no more wrong than the current O2 in
   every world: intact,
   paralog-annotated isoforms, an uncatalogued in-genome copy, a removed copy, silent copies. It recovers >= half of the
   uncatalogued copies' reads, drops >= 80% of silent copies and <= 10% of expressed ones, and shrinking the universe does not
   cost correct calls." Any failing bar: not claimed.
2. **J6-J7 pass:** claimed, "integrated O3 detects >= 75% of detectable removed copies with at most one false call in 32
   intact, paralog-annotation, uncatalogued-copy and silent-copy families". J7 not measurable: claimed without the
   removed-copy clause.
3. Each claim stands or falls as a whole.

Scope: simulated reads (curated CAT isoforms, HiFi errors, heterozygous sites at 0.1%, uniform expression) on the real CHM13
genome. Real data's other species, allelic imbalance and expression variation are not modelled.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (held out, 2026-10-16)

Run: code and prereg `1e513942`, the frozen copy_assign binary, seed 81, `--het 0.001`. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout7/` (`integ.json`, `o3_score.json`, `o3_summary.txt`); summaries from `summary.py integ`
and `summary.py o3`.

**O2** (correct / wrong; unknowable placements in brackets):

| world | current O2 | padded union O2 | **integrated O2** | ratio of correct calls to the current O2's |
|---|---|---|---|---|
| INTACT | 2,360 / 0 | 2,476 / 0 | **2,476 / 0** | x1.049 |
| CROSS | 2,697 / 0 | 2,860 / 15 | **2,955 / 11** | x1.096 |
| OUTSIDE | 2,081 / 68 | 2,193 / 193 | **2,533 / 180** | x1.217 |
| REMOVED | 2,081 / 45 (n/a) | 2,193 / 0 (193) | **2,193 / 0 (193)** | x1.054 |
| SILENT | 1,723 / 0 | 1,762 / 0 | **1,834 / 0** | x1.064 |

- **OUTSIDE.** O1 seeded 7 new copies, all on a truth copy, and 348 of the 481 uncatalogued reads went to their own copy. The
  integrated O2's 180 wrong calls are:
  - 60 + 60: two uncatalogued copies that were not seeded (CTAGE.2, DDX11L.9), whose reads went to a near-identical paralog;
  - 60: WASH.7's reads called to WASH.7 under the DDX11L family. The WASH and DDX11L families are adjacent subtelomeric genes,
    and the uncatalogued WASH.7 was seeded in both, since a new locus is off every family's catalog.
- **CROSS.** 11 paralog-structure reads were called to a wrong copy: DDX11L 8, WASH 3.
- **SILENT.** O1 dropped 10 of the 14 silent copies and 0 of 220 expressed copies over the five worlds. The 4 kept silent
  copies had >= 5 primary alignments: mapq-0 coin flips of identical expressed copies. Keeping them is the rule by design (an
  identical expressed copy must never be dropped).
- **REMOVED.** The removed copies' median `de`: TUBB8 0.087, SPRR2 0.076, WASH 0.0124 (detectable); AQP7 0.0077, FAM74A
  0.0093, FAM138 0.0043, CTAGE 0.0023, DDX11L 0.0022 (unknowable).

**O3** (8 families per world):

| | CUR | S23 | H | **I** |
|---|---|---|---|---|
| REMOVED detected (3 detectable) | 1 (TUBB8) | 2 | 2 | **2** (SPRR2, TUBB8) |
| false-call families, INTACT + CROSS + OUTSIDE + SILENT (of 32) | 8 | 0 | 0 | **0** |

WASH (median `de` 0.0124, near delta) is missed by every arm.

| bar | value | |
|---|---|---|
| J1 | lowest world ratio 1.049 | pass |
| J2 | OUTSIDE 180 > 68; CROSS 11 > 0 | **fail** |
| J3 | 348 / 481 = 72% | pass |
| J4 | silent copies dropped 10 / 14 = 71%; expressed 0 / 220 | **fail** |
| J5 | 1,834 vs 1,762 | pass |
| J6 | 0 of 32 | pass |
| J7 | 2 of 3 detectable (0.67), CUR's one included | **fail** |

#### Verdict (2026-10-16)

**O2: not claimed (J2 and J4 fail).** On this fresh set the integrated O1 -> O2:
- assigned 5-22% more reads correctly than the current O2 in every world (J1), with no wrong call in INTACT, REMOVED and
  SILENT;
- recovered 72% of the uncatalogued copies' reads (J3);
- shrank the universe without losing correct calls (J5);
- dropped no expressed copy.

But:
- **wrong calls:** it made more than the current O2 where uncatalogued copies are near-identical to catalogued ones and
  were not seeded (CTAGE, DDX11L), and where an uncatalogued copy was seeded into two adjacent families (WASH / DDX11L);
- **silent copies:** identical silent copies keep coin-flip support, so only 71% were dropped.

**O3: not claimed (J7 fails, 2 of 3).** No false call over 32 families. WASH's removed copy, near delta, was missed by
every arm, as was round 3's case near the boundary.

**What the integration showed, across development and this held-out set:**
- **What it fixes:** RNA-supported seeding fixes O2's main failure mode where it can, namely uncatalogued copies with
  uniquely placed reads. On this set wrong calls went from 193 for padded union to 180, with 348 reads now right.
- **What it does not fix:** copies near-identical at mapq 0 (unseedable), and alleles.
- **A bug it exposed:** a seeded locus must belong to one family only. Fixing that is a method change for a new prereg.

**Development O3 with the final seeding (promised in "Known before running").** Arm I made no false call in any of the 12
development worlds. It detected every detectable removed copy: round 6's set 3 of 3 (CUR 1), round 4's set 6 of 6 (CUR 5).

**End to end, post hoc (`summary.py e2e`, `read_fate`).** Every read's fate under each system, by the kind of copy it came
from, with an absent-copy call taking precedence over an assignment. On this held-out set:

| world | current: correct / error / abstain | integrated: correct / error / abstain |
|---|---|---|
| INTACT | 2,348 / 44 / 615 | 2,476 / 0 / 531 |
| CROSS | 2,682 / 43 / 943 | 2,955 / 11 / 702 |
| OUTSIDE (catalogued + uncatalogued reads) | 2,069 / 233 / 705 | 2,533 / 120 / 354 |
| REMOVED (+ unknowable placements) | 2,190 / 44 / 473 (+68) | 2,405 / 13 / 325 (+193) |
| SILENT | 1,705 / 46 / 415 | 1,834 / 0 / 332 |

- The current method's errors are mostly its absent-copy calls on present copies' reads, which the O2 bars did not count.
- Reads of the detectable removed copies end inside an absent-copy call: 181 / 181 for the integrated pipeline, 121 for the
  current method.
- On the 12 development worlds the integrated pipeline has more correct fates and no more errors in 9.
  - Round 4's INTACT and CROSS: 5 and 12 fewer correct reads, with 0 errors either way.
  - Round 6's REMOVED: 227 more correct reads and 1 more error.

This is computed after the results were seen, so it is not a claim. It is the basis for a pre-registered end-to-end test on a
new family set.


---

## PREREG_o3_graph_round8_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round8_2026-10-16.md`)*

## Round 8: the integrated O1 -> O2 -> O3, end to end, on a ninth handpicked family set (2026-10-16)

Pre-registration, written before any world of the ninth set was built: no read of these families has been simulated. Code:
branch `o3-e2e` (`bench/o3_graph/`), frozen at the commit that adds this file.

**User request (2026-10-16):** "fix the bug so an added copy belongs to only one family. Then pre-register an end-to-end test
on a fresh 9th family set, with bars on read outcomes per world type (correct ≥ current, errors ≤ current)".

**Question.** On a fresh family set, does the integrated pipeline do at least as well as the current method in every world
type? That means at least as many reads with a correct outcome (on their own copy, or inside a call for their absent copy) and
no more reads with a wrong one.

### The two systems

**Current** (unchanged since round 1):
- copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`) on the curated catalog;
- its `--psv-clusters` O3: a cluster with status `unexplained` and estimated divergence > delta (0.00958) is an absent-copy
  call.

**Integrated** (round 7 plus the one-family rule):
- **O1**, `o1seed.py`:
  - keep each curated copy on which >= 5 family reads have their primary alignment (+-1 kb), and drop the others (silent);
  - add off-catalog loci with >= 5 uniquely placed reads;
  - give each copy the intron chains of near-unique reads (mapq >= 50) that it owns;
  - **new (1e3960fe):** a new locus found in several families stays only in the family whose curated copies its reads align
    to most (aligned bases, unpadded copy spans). A tie seeds it in none.
- **O2**, `o2graph.py run W integ 2`: padded union O2 on the seeded catalog.
- **O3**, `detect.py` arm I: the pooled clusters of reads that O2 leaves unexplained, with the mosaic rule, branch calls and
  the genome check.

Both run per world from one driver, `bench/o3_graph/round8.sh`.

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout8`, 1d391d5d):
- NOTCH2NL, SRGAP2, GH (GH1/2, CSH1/2, CSHL1), USP17L, TSPY, RBMY, BPY2, CDY, OR4F, LCE1;
- KRTAP5, PNMA6, PPIAL4, FCGR1, NOMO, CROCC, RGPD, FAM86, GPR89, POM121.

None of them was in any earlier world. FAM86, GPR89 and POM121 were round-4 candidates past its first 8.

The build rule is the same as in rounds 1-7:
- copies are linked at identity >= 0.90 over >= half the shorter canonical transcript;
- genes span <= 60 kb;
- a seeded sample of 10 copies is taken when a family has more.

NOTCH2NL and SRGAP2 form no family under that rule. The set is the first 8 with >= 3 copies:

| family | genes | copies (of linked) / isoforms |
|---|---|---|
| GH | CSH1, CSH2, CSHL1, GH2 (chr17) | 4 / 16 |
| USP17L | USP17L3, L4, L23, L25-27 (chr4, chr8) | 10 (of 59) / 10 |
| TSPY | TSPY1-4, 8, 10 (chrY) | 10 (of 35) / 38 |
| RBMY | RBMY1A1, 1D, 1E, 1F (chrY) | 5 / 16 |
| BPY2 | BPY2, BPY2B (chrY) | 3 / 12 |
| CDY | CDY1, CDY2A, CDY2B (chrY) | 4 / 4 |
| OR4F | OR4F4, OR4F5 (chr3, chr11, chr15) | 3 (of 8) / 3 |
| LCE1 | LCE1D, E, F (chr1) | 3 (of 6) / 5 |

Files: `/mnt/linuxdisk/tmp/o3_graph/families8.json` (round 7's file plus this set), `families8_all.json` and `families8.log`.
The rebuild reproduces every earlier set exactly.

**This is the hardest set so far.** Six of the eight families have near-identical copies:
- the four Y-ampliconic families;
- OR4F;
- the USP17L tandem array.

Their reads tie at mapq 0, the failure mode that no step of either system fixes.

### Worlds

All with heterozygous alleles (`--het 0.001`, as in round 7), seed **93**:
- INTACT;
- CROSS (`--projection liftover`);
- OUTSIDE;
- REMOVED (the special copy hard-masked, with its own masked genome and index);
- BRANCH (an extra absent copy with a novel 150-bp exon);
- SILENT.

**Seed (detectability pre-check, as in round 7).** `precheck.py families8.json heldout8 91 101` (`precheck8.txt`) gives seed 93,
the smallest seed in 91-100 with >= 3 detectable removed copies:
- detectable: GH 0.0367, USP17L 0.0163, LCE1 0.0391;
- within delta: TSPY 0, RBMY 0.0005, BPY2 0, CDY 0.0004, OR4F 0.

**Commands:**
```
R=/mnt/linuxdisk/tmp/o3_graph/real/heldout8; F=/mnt/linuxdisk/tmp/o3_graph/families8.json
python3 bench/o3_graph/realworld.py mask $R/masked $F --set heldout8 --seed 93
minimap2 -x splice:hq -t 4 -d $R/masked/masked.splice.mmi $R/masked/masked.fa
for sc in INTACT CROSS OUTSIDE REMOVED BRANCH SILENT; do bench/o3_graph/round8.sh $R $F heldout8 93 $R/masked $sc; done
python3 bench/o3_graph/summary.py e2e $R/e2e.json $R/INTACT $R/CROSS $R/OUTSIDE $R/REMOVED $R/BRANCH $R/SILENT
python3 bench/o3_graph/summary.py bars $R/e2e.json
```

### Read outcomes (frozen: `summary.py` `read_fate`, `cmd_e2e`)

**Class.** Each read is classed by the copy it came from:
- catalogued;
- uncatalogued but in the genome (OUTSIDE);
- absent and detectable: BRANCH's extra copy, or a REMOVED copy whose reads' median `de` on the masked genome is > delta;
- absent and undetectable.

**Fate.** Under each system, a read's fate is one of:
- **correct:**
  - inside an absent-copy call, for an absent copy's read;
  - otherwise, assigned to its own copy (an uncatalogued copy can only be its own once O1 has seeded it).
- **error:**
  - inside an absent-copy call, for a present copy's read;
  - assigned to another copy;
  - a detectable absent copy's read assigned to a present copy.
- **unknowable:** an undetectable absent copy's read assigned to a present copy.
- **abstain:** no call.

An absent-copy call takes precedence over an assignment.

**A read called in several families** is judged by its call in its own family, as in round 7 and the overview. This favours
the current method: over the 43 earlier worlds, 115 reads had calls in two families, all of them from the current method.

### Bars, per world, all six worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Correct and errors are summed over every read class. Abstain and unknowable count for neither. `summary.py bars` computes both
bars (`round8_bars`, tested in `test_summary.py`).

### Reading

- **PASS** (B1 and B2 hold in all six worlds). Claim: "on a fresh family set, the integrated O1 -> O2 -> O3 pipeline gives,
  in every world type, at least as many correct read outcomes as the current method and no more wrong ones".
- **FAIL** otherwise. No partial claim. The failing worlds and their margins are reported.

### Known before running

- **The 43 earlier worlds** (rounds 1-7, post hoc, `docs/O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md`): both bars held in 38.
  - B1 failed by 1, 2, 5 and 12 reads in round 3's and round 4's INTACT and CROSS, where neither system made an error.
  - B2 failed by 1 error in round 6's REMOVED, alongside 227 more correct reads.
  - A strict bar can therefore fail by a handful of reads.
- **Round 7's held-out set** passes these bars post hoc, although its pre-registered per-objective O2 bars failed.
- **The one-family rule on development worlds.** These are the six worlds where a new locus had been seeded twice: set 1's
  five worlds (NPIPA/NPIPB) and round 7's OUTSIDE (WASH/DDX11L).
  - It was never worse.
  - Set 1's OUTSIDE errors fell from 30 to 17, and set 1's CROSS gained 10 correct reads.
  - Nothing else changed.
- **BRANCH with alleles** had never been run before. On round 7's families (development):
  - current: 2,406 correct and 134 errors;
  - integrated: 2,888 correct and 9 errors;
  - arm I made 16 calls, all true.
- **The integrated pipeline's false absent-copy calls** over the 43 worlds were all one cluster (NPIPA.5's reads hosted on
  NPIPB, at delta). The fix does not change it.
- **Its errors on uncatalogued copies** (291) come from near-identical copies whose reads tie at mapq 0. This set has six such
  families.

### Reported beside the bars (not gating)

- the fates per read class, and the error sources (as in the overview);
- **O1:** copies kept and dropped (silent versus expressed), new copies (on a truth copy or not), and loci resolved by the
  one-family rule;
- **O3:** arm I's calls, true or false;
- the same per family.

### After the run

- An independent review by a fresh subagent: the code path, the numbers, and this document's reading.
- Register row 1306.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-16; worlds built 13:22-13:31 with the code of e4d0703b)

### Verdict: FAIL. B2 fails in CROSS by one read. Not claimed.

`summary.py bars heldout8/e2e.json` (frozen scorer). An independent review reproduced every cell with its own scorer, written
from this document's text.

| world | correct: current / integrated | B1 | errors: current / integrated | B2 |
|---|---|---|---|---|
| INTACT | 1,085 / 1,293 (x1.19) | yes | 0 / 0 | yes |
| CROSS | 1,226 / 1,427 (x1.16) | yes | 0 / **1** | **no** |
| OUTSIDE | 964 / 1,222 (x1.27) | yes | 301 / 60 | yes |
| REMOVED | 1,144 / 1,345 (x1.18) | yes | 1 / 0 | yes |
| BRANCH | 1,119 / 1,690 (x1.51) | yes | 95 / 9 | yes |
| SILENT | 724 / 932 (x1.29) | yes | 0 / 0 | yes |

The six worlds share most of their reads, so their totals are not independent counts:
- INTACT, OUTSIDE and REMOVED simulate the same reads;
- CROSS and BRANCH contain all of INTACT's reads;
- SILENT's reads are a subset of them.

Each world is the unit of comparison.

**The one error** is read s0001190, from TSPY.2 on haplotype 1, which has no alleles. Its isoform is projected from TSPY.7's
annotation (CROSS). The integrated O2 called it TSPY.1. The review replayed the simulator and found the cause:
- One of the read's two simulated substitutions (read position 3) hits the only base that separates TSPY.1 from TSPY.2 over
  its 663-bp structure, and turns it into TSPY.1's base.
- minimap2 soft-clips the first 4 bases on TSPY.2, and O2's distance (`_paf_dist`) counts clipped bases as edits. That gives
  1 edit against 5, which is past the margin of 2.
- With end-to-end alignment, the distance to TSPY.2 would be 2 and the read would abstain.

The current method leaned the same way (copy 1, one decisive site, p 0.001) but abstained as ambiguous. The result stands as
run; the frozen code is what was tested. A semi-global O2 distance is a candidate change, to be tested on a fresh set.

### Beside the bars (not gating; corrected by the review)

**Where the errors are:**
- current:
  - 301 reads of uncatalogued copies. 180 of them are its absent-copy calls on copies that are in the genome; it has no
    genome check. Without those it still has 121, against 60;
  - 95 reads of novel-exon absent copies, placed on a present copy;
  - 1 read of a removed copy;
  - no error on a catalogued copy's read.
- integrated:
  - 60 reads of uncatalogued copies;
  - 9 BRANCH reads, all TSPY.B reads without the novel exon, which are essentially the parent copy's sequence;
  - the CROSS read.

**O1:**
- dropped 7 of 12 silent copies, and none of the expressed copies in any world;
- seeded 3 of the 8 uncatalogued copies, but only GH.3 and LCE1.2 received their own reads:
  - RBMY.0's seeded copy was built from 55 of its reads, but its observed structures miss the reads' starts;
  - its 61 reads end `thin` (52) or `tied` (9);
- the other five uncatalogued copies (BPY2.1, CDY.0, OR4F.2, TSPY.3, USP17L.8) are near-identical to catalogued copies, and
  their reads tie;
- in REMOVED it seeded OR4F17 (chr19), which is none of the family's copies. The masked OR4F.2's reads (undetectable) land
  there uniquely;
- **the one-family rule resolved no locus**, so round 8 is no held-out evidence for it. Its evidence remains the development
  worlds.

**O3 (arm I):**
- 20 calls, none false;
- REMOVED: all 181 reads of the three detectable removed copies (GH.3, LCE1.2, USP17L.8) are inside pure calls. Under the
  single-call rule of rounds 4-7 (one call holding >= 80% of the copy's reads), 2 of 3 are detected: GH.3's 61 reads are split
  over five calls, and the largest holds 43 (0.70);
- BRANCH: **7** novel-exon copies were built, not 8. USP17L's special copy has a single-exon main isoform, so it has no
  middle to insert an exon into. All 7 have a call made mostly of their reads. Under `score.py`'s purity rule (0.9) it is 6
  of 7: OR4F.B's only call holds 34 novel-exon reads plus 26 of the same copy's reads without the exon (0.57).

### Scorer defect found by the review (bars unchanged)

`cmd_e2e` scored a call on a seeded copy that matches no truth copy as an abstention, not as a placement. In REMOVED, 59 of
OR4F.2's reads placed on OR4F17 should therefore be `unknowable`, not `abstain`.

The fix (`integrated_calls`, with a test) counts such a call as a placement on a locus that is no copy. Re-scoring all 50
worlds with it changes only REMOVED worlds, and only abstain to unknowable:
- set 2: 28 reads;
- set 6: 48 reads;
- this set: 59 reads.

No correct or error count changes, so neither bar changes. The verdict above uses the frozen scorer.

### Corrections to the text above

The family table's "(of 8)" for OR4F and "(of 6)" for LCE1 count candidate genes, not linked copies: their linked components
have 3 copies each.

### Reading

**The test does not support the pre-registered claim**, and no partial claim is made.

It fails on one read: an integrated O2 call among near-identical copies. A single sequencing error at the one site that
separates TSPY.1 from TSPY.2, plus counting soft-clipped bases as edits, pushed that call past the margin.

Everything else matches the 43 earlier worlds:
- the integrated pipeline gets more reads right in every world type;
- it has fewer errors in five of six;
- its O3 makes no false call.

A per-world strict error bar is decided by single reads wherever the current method makes no error (here INTACT, CROSS and
SILENT). The current method makes none there because it abstains more.


---

## PREREG_o3_graph_round9_2026-10-16.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round9_2026-10-16.md`)*

## Round 9: the integrated pipeline v2 under realistic reads, end to end, on a tenth set of full families (2026-10-16)

Pre-registration, written before any world of the tenth set was built: no read of these families has been simulated. Code:
branch `o3-fast` (`bench/o3_graph/`), frozen at the commit that adds this file.

**Why (the user's goal, 2026-10-16):** "integrate O1, O2 and O3 more so it's faster and cheaper to test full families; also
improve the simulations to account for more realistic ones and determine if PSVs can also include SUNs and probably some
influence of the cigar string, nm, de tags".

The development work behind this test is `docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` (sets 8 and 9):
- **Calibrated realistic reads remove the integrated pipeline's advantage on catalogued copies.**
- **CIGAR-aware evidence restores it:**
  - O2 counts indel EVENTS, not indel bases, and aligns clipped ends against the path;
  - O3 drops, read by read, the members the genome explains.

**Question.** On a fresh set of full (uncapped) families, with reads from two calibrated real libraries, does the integrated
pipeline v2 do at least as well as the current method in every world type? That means at least as many correct read outcomes
and no more wrong ones.

### The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v2:**
- **O1:** `o1seed.py` (round 8's, unchanged);
- **O2:** `o2graph.py run W integ 2 --distance gapcomp`;
- **O3:** `detect.py W --member-genome-check`, arm I.

Both run per world from `fastset.py run ... --o2-distance gapcomp --member-genome-check`, with the set's local reference
(`localref.py`, validated in the development doc).

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout9`):
- PPIAL4, FCGR1, PGA, H4C, LCE2, LCE3, KRTAP4, AGAP;
- ELOA3, FAM25, CYP3A, FAM47, FAM9, KRTAP9, DAZ.

PPIAL4 and FCGR1 were round-8 candidates past its first 8, never built into a world; none of the others was in any earlier
list.

The build rule is rounds 1-8's (identity >= 0.90 over >= half the shorter canonical transcript, span <= 60 kb), with
**no copy cap** (`--cap 1000`: full families). The set is the first 8 with >= 3 copies (H4C and LCE3 link 2 copies each):

| family | copies | isoforms | chromosome |
|---|---|---|---|
| PPIAL4 | 7 | 7 | chr1 |
| FCGR1 | 3 | 9 | chr1 |
| PGA | 3 | 12 | chr11 |
| LCE2 | 3 | 3 | chr1 |
| KRTAP4 | 5 | 5 | chr17 |
| AGAP | 7 | 15 | chr10 |
| ELOA3 | 12 | 12 | chr18 |
| FAM25 | 3 | 4 | chr10 |

Files: `/mnt/linuxdisk/tmp/o3_graph/families9.json` (round 8's sets plus this one), `families9_all.json`, `families9.log`.

### Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under two read profiles:
`bench/o3_graph/profiles/A119b_testis.json` and `bench/o3_graph/profiles/public_testis.json`. That makes 12 worlds.
- **Profiles:** calibrated truncation, homopolymer-biased indels, the donor's het and hom alleles, retained introns, and
  per-copy expression spread. With a profile, `--het 0` means the donor's alleles come from it.
- **Seed 101.** `precheck.py families9.json heldout9 101 111` (`precheck9.txt`) gives the smallest seed in 101-110 with >= 3
  detectable removed copies:
  - detectable: LCE2 0.068, KRTAP4 0.039, FAM25 0.017, AGAP 0.014;
  - within delta: PPIAL4, FCGR1, PGA, ELOA3.

**Commands** (one per profile; the local reference is built by the first run):
```
F=/mnt/linuxdisk/tmp/o3_graph/families9.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout9
V="--o2-distance gapcomp --member-genome-check --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout9 101 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout9 101 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis.json
python3 bench/o3_graph/summary.py e2e $O/<profile>/e2e_pr.json <its six worlds> --per-read-detectability
python3 bench/o3_graph/summary.py bars $O/<profile>/e2e_pr.json
```

### Read outcomes

These are round 8's (`summary.py` `read_fate`, `cmd_e2e`), with one change, made on development data before this
registration: **a REMOVED copy's read is detectable by its own primary `de` on the masked genome (> delta)**
(`--per-read-detectability`).
- **Why:** on set 9 the one remaining error under both profiles was a read of a detectable copy (median `de` 0.018) that is
  near-identical to a present copy over its own truncated span. No method can tell that read apart.
- **The same rule applies to both systems.**
- **Secondary report:** round 8's copy-level classification (`e2e.json`) is reported too, not gating.

A read called in several families is judged by its call in its own family, as before.

### Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither (`summary.py bars`).

### Reading

- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "with reads calibrated on two real libraries, the integrated O1 -> O2 ->
  O3 pipeline v2 gives, in every world type of a fresh set of full families, at least as many correct read outcomes as the
  current method and no more wrong ones".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.

### Known before running (development sets 8 and 9)

- **Integrated v1** (round 8's code) fails under both profiles on both sets: see the development doc, section 3.
- **Integrated v2, set 9:**
  - B1 holds in all 12 worlds;
  - copy-level scoring: B2 fails only in REMOVED, 1 vs 0, under both profiles;
  - per-read detectability: PASS under both profiles.
- **Integrated v2, set 8:**
  - public profile: PASS under both scorers;
  - A119b profile: FAIL under both. B1 fails in INTACT, 2,215 vs 2,234 correct. Under copy-level scoring B2 also fails in
    REMOVED, 5 vs 3; per-read scoring makes it 3 vs 6.
  - The INTACT shortfall is O2's `unexplained` reads: 37 in SPRR2, 12 in AQP7. The current method assigns them.
  - A MAPQ-60 genome fallback for them was measured and rejected: 17 wrong of 329 outside BRANCH, 159 wrong of 220 in BRANCH.
- **Development summary:** under per-read scoring v2 passes 3 of the 4 realistic set x profile runs (all but set 8, A119b);
  v1 passed none.
- **Integrated v2 under idealized reads** (copy-level scoring): PASS on sets 8 and 9. On set 9 it fixes round 8's single
  CROSS error, at a cost of 4-6 correct reads per world. On set 8 it gets fewer reads right than v1 (INTACT 2,407 vs 2,476)
  with fewer errors (CROSS 3 vs 11, REMOVED 0 vs 13).
- **What was chosen on these sets:** the options (`gapcomp`, `--member-genome-check`) and the per-read rule. This tenth set
  is their first test.
- **A per-world error bar is decided by single reads** wherever the current method makes no error.

### Reported beside the bars (not gating)

- the fates per read class, and the error sources;
- O1's drops and new copies;
- arm I's calls, true or false;
- **the run's cost:** `timing.json`, every step's wall time and peak RSS.

### Amendment 1 (2026-10-16 15:28, before any result)

**The first run crashed in O3 under both profiles.**
- `vg msga` (v1.73, a deprecated subcommand) fails an assertion (`find_mems_deep: lcp_maxima.size() == mems.size()`) on
  KRTAP4's five loci in FASTA order, both by default and with `-a`.
- Any two of the loci build, and every other order tried builds.

**The fix (`detect.py`, `msga_orders`, tested):**
- retry `vg msga -a` with the loci reversed, then each rotation;
- if every order fails, the family gets no graph and no branch calls; its PSV clustering still runs.
- The reversed order builds KRTAP4's graph.

**Nothing else changes.**
- On development worlds, where msga never failed twice, `detect.py`'s output is byte-identical.
- No outcome of round 9 was computed before this amendment: both runs stopped at O3, and no e2e or bars were written.
- The two crashed run directories were deleted; the run is repeated from the start with the same commands.

### After the run

- An independent review by a fresh subagent.
- Register row 1307.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-16; run after amendment 1 with the code of 8235f733)

### Verdict: FAIL under both profiles. Not claimed.

Bars from `summary.py bars e2e_pr.json` (per-read detectability). Copy-level scoring gives the same verdict, cell for cell
except REMOVED's error counts (A119b 10 vs 2, public 5 vs 0).

| world | A119b profile: correct (current / integrated) | errors | B1, B2 | public profile: correct | errors | B1, B2 |
|---|---|---|---|---|---|---|
| INTACT | 1,149 / 1,109 | 2 / 0 | **no**, yes | 994 / 976 | 1 / **5** | **no**, **no** |
| CROSS | 1,264 / 1,227 | 2 / 1 | **no**, yes | 1,085 / 1,071 | 2 / **5** | **no**, **no** |
| OUTSIDE | 868 / 1,069 | 177 / 3 | yes, yes | 740 / 944 | 174 / 8 | yes, yes |
| REMOVED | 1,038 / 1,013 | 11 / 3 | **no**, yes | 917 / 964 | 4 / 0 | yes, yes |
| BRANCH | 1,181 / 1,301 | 219 / 55 | yes, yes | 994 / 1,108 | 180 / 67 | yes, yes |
| SILENT | 698 / 686 | 2 / 0 | **no**, yes | 552 / 542 | 1 / **2** | **no**, **no** |

**Where the integrated pipeline wins.** On a fresh set of full families with realistic reads, it is far ahead wherever the
catalog is incomplete:
- OUTSIDE: +201 and +204 correct, 3 and 8 errors against 177 and 174;
- BRANCH: +120 and +114 correct, 55 and 67 errors against 219 and 180;
- REMOVED (public): +47 correct, 0 errors against 4.

**Where it loses: catalogued copies.** It gets 1-3.5% fewer reads right in INTACT, CROSS and SILENT, and in REMOVED under
A119b.
- **Fewer correct reads:** the losses are reads O2 leaves `unexplained`, spread over the families. For INTACT under A119b:
  AGAP 15, KRTAP4 12, LCE2 11, FCGR1 9, PPIAL4 5. The current method assigns them.
- **The public profile's extra errors:**
  - INTACT's 5 are all PGA reads (pepsinogen A3/A4/A5, near-identical copies) called to the wrong copy by O2;
  - CROSS and SILENT were not broken down.
- This is the same mechanism development set 8 showed under A119b (INTACT -19).

**Cost.** 168 s at 0.76 GB peak (A119b) and 149 s at 0.74 GB (public) for six worlds of 43 copies each. The local reference
was built once, in the first run.

### Reading

**The pre-registered claim fails under both profiles.** No partial claim is made.

Three determinations hold on the held-out set as on development data:
- **Realistic reads change the comparison:** the advantage that idealized reads showed on catalogued copies does not survive.
- **The integrated pipeline's strength is the incomplete catalog:** copies the catalog lacks (OUTSIDE), copies absent from the
  reference with a novel exon (BRANCH), and their errors.
- **The remaining weakness is O2's divergence gate on realistic reads.** These are reads with retained introns, long indels or
  short spans that no annotated path explains. A genome-placement fallback for them was measured on development data and
  rejected for its errors.

### Independent review (2026-10-16) and corrections

**The verdict holds.**
- The reviewer's own scorer, written from this document's text, reproduces `e2e.json` and `e2e_pr.json` exactly for both
  profiles.
- The committed code rebuilds all 12 worlds byte for byte. Re-running O1, O2 (`gapcomp`) and O3 (`--member-genome-check`) on
  copies reproduces their outputs.
- The families share no gene with any earlier set. The seed came from the pre-check before any run. The final runs started
  five seconds after the amendment commit.

**Corrections to the results above:**

1. **The public profile's extra errors are an O1 step, not an O2 one.**
   - All 12 integrated errors under the public profile (INTACT 5, CROSS 5, SILENT 2) are reads of PGA.2 (PGA5).
   - O1 dropped PGA.2: it had 5 reads but only 4 primaries on its locus, 2 of them uniquely placed, so it fell under the 5-read
     support rule.
   - O1 dropped 2-4 expressed copies with 4-5 reads each in every public-profile world (AGAP.0, LCE2.1, PGA.2, PPIAL4.3), and
     none under A119b.
   - Those drops account for 10 of INTACT's 18-read deficit, 6 of CROSS's 14, and all 10 of SILENT's.
   - The verdict still fails through B1 in INTACT and CROSS, and through every A119b failure.
   - Cause: the public profile's expression spread (1.16) was inflated by zero-read copies and Poisson counting noise (item 3).
     It produced copies with 4-5 reads, and O1's fixed support rule drops them.
2. **The pre-registered side reports:**
   - **Arm I:** no false call in any of the 12 worlds.
   - **REMOVED:** under A119b arm I called 2 of the 4 detectable copies (LCE2.0, KRTAP4.0, plus PPIAL4.4). It missed AGAP.5 and
     FAM25.2, whose median `de` is 0.0149 and 0.0164, although O2 left 54 of 57 and 21 of 21 of their reads unexplained.
     Under the public profile it called all four.
   - **BRANCH:** 6 of 6 novel-exon copies under A119b, 5 of 6 under the public profile.
   - **O1:** dropped 8 of 13 silent copies under both profiles; the expressed-copy drops are in item 1.
   - **"REMOVED losses are reads O2 leaves unexplained" is not the whole story.** The current method caught FAM25.2, which
     arm I missed; that accounts for 15 of the 55 correct reads lost in A119b's REMOVED.
3. **The expression spread.** Re-estimated over copies with reads, less Poisson noise
   (`calibrate.py` `expression_log_sd`): A119b 0.86 (unchanged), public 1.03 (not 1.16).
4. **Cost.** The per-run peaks (0.74-0.76 GB) leave out the one-time local reference build. That build loads the whole-genome
   index, 38-75 s at about 14 GB, once per family set. A new family set therefore still peaks near 14 GB, against 19.4 GB
   before. Round 9's own build ran inside the deleted crashed run.
5. **"Realistic reads change the comparison" was not tested on this set:** there was no idealized-read run of set 10. A post
   hoc idealized run is reported in the development doc.
6. **Smaller points:**
   - The local reference was validated against whole-genome mapping on idealized reads only. minimap2's repeat filter differs
     on the local reference (99 vs 969 occurrences), so realistic reads are untested there.
   - minimap2 breaks MAPQ-0 ties partly on the read NAME, so mapping by sequence moves tie choices; the bars were identical.
   - Before amendment 1, the crashed runs' current, O1 and O2 outputs existed and were deleted. The order retry cannot be tuned
     to the bars.
   - Per-read detectability partly reflects sequencing noise.
   - A BRANCH read without the novel exon still counts as detectable.


---

## PREREG_o3_graph_round13_2026-10-17.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round13_2026-10-17.md`)*

## Round 13: the integrated pipeline v5', independent of the current method, on a fourteenth set of full families (2026-10-17)

Pre-registration, written before any world of the fourteenth set was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal: integrate O1, O2 and O3 so that full families are fast and cheap to test, with realistic
simulations. Round 12 (register 1310) failed under both profiles. Its review found two problems:
- **Most of its error reads are structural.** Retained introns are absorbed as one gap-compressed edit, and two copies are
  DNA-identical.
- **O3's cost was set by the branch arm,** whose all-against-all alignment ran away when GraphAligner failed on the 80-copy
  DUX4 graph.

v5' removes the cost problem. It does not remove round 12's O2 errors (see "Known before running").

**Question.** On a fresh set of full families, including a 43-copy pseudogene array and retrocopy families, does v5' do at least
as well as the current method in every world type, with both calibrated v2 profiles? That means at least as many correct read
outcomes and no more wrong ones.

### The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v5'** = round 12's v5 + `--branch-genome-novel` (16ba2508):
- **O1:** `o1seed.py W --twins`.
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `run W integl 2 --error-floor`.
- **O3:** `detect.py W --member-genome-check --member-attribution --error-floor --o1-hosts --branch-genome-novel --arms I`.
  - It is hosted by O1 and reads nothing of the current method's run.
  - The branch arm clusters only reads whose own genome alignment leaves >= `MIN_NOVEL` bases unexplained. The genome check
    applies the same criterion to every call afterwards. A read the genome explains carries no novel exon, and this bounds
    the all-against-all when GraphAligner fails on a graph.
- **Measured and not adopted:**
  - **`--distance gapstruct`** (gap compression only below 30 bp). It still fails round 12, adds an error in set 8 A119b
    INTACT and SILENT, and costs set 10 A119b REMOVED 3 correct reads.
  - **`--mapped-sites`.**

```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --heavy-current
    --step-timeout 1800 --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```
`--heavy-current` and `--step-timeout 1800` are set from the start: round 12's 80-copy array needed both, and outputs depend on
neither.

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout13`, 5c1f3ab9):
- **MRGPRX:** a round-12 candidate with >= 3 copies past its first 8, never built into a world.
- **S100A7, KRTAP12, KRTAP13, KRTAP20, KRTAP21, KRTAP6, HLA-ABC, H1, OR2A, SEPTIN7P, CTSLP, UBE2Q2P, PGAM1P, OR52, OR7E:** in
  no earlier list.
- **Spares, appended before any run (bb60c046)**, because the first build gave only 5 families with >= 3 copies: YWHAZP,
  PTMAP, TPT1P, FTH1P, EEF1A1P, KRT8P (retrocopy families, smallest first).

**Build rule:** round 12's, `families.py build families13_all.json --cap 1000 --link locus --min-cover 0.8`, every group. The set
is the first 8 families with >= 3 copies, in list order:

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| MRGPRX | 4 | 5 | chr11 |
| SEPTIN7P | 5 | 8 | chr7, chr9 |
| CTSLP | 6 | 8 | chr9, chr10 |
| UBE2Q2P | 4 | 8 | chr15 |
| OR7E | 43 | 44 | chr2, chr3, chr4, chr7, chr8, chr9, chr10, chr11, chr13, chr14 |
| YWHAZP | 5 | 5 | chr2, chr6, chr9, chr10 |
| PTMAP | 9 | 9 | chr3, chr5, chr7, chr9, chr12, chr13, chr20 |
| FTH1P | 13 | 13 | chr1, chr2, chr3, chr5, chr6, chr8, chr9, chr11, chr13, chrX |

89 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families13.json` (this set only), `families13_all.json`, `families13.log`.

### Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under the two v2 profiles,
`profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json` (`--het 0`). That makes 12 worlds.

**The seed:** `precheck.py families13.json heldout13 181 191` (`precheck13.txt`; 70.5 s, 14.5 GB) gives the smallest seed in
181-190 with >= 3 detectable removed copies: **seed 181**.
- Detectable: MRGPRX 0.0809, YWHAZP 0.0597, OR7E 0.051, SEPTIN7P 0.0436, FTH1P 0.0435, UBE2Q2P 0.0433, CTSLP 0.0102.
- Within delta: PTMAP 0.0061.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families13.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout13
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout13 181 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout13 181 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```

### Read outcomes, bars and reading

These are round 12's, unchanged:
- **Read outcomes:** per-read detectability less the world's error floor. Round 10's rule and copy level are secondary.
- **Bars:** B1 (integrated correct >= current) and B2 (integrated errors <= current), per world, all 12 worlds.
- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "the integrated O1 -> O2 -> O3 pipeline v5', which uses nothing from the
  current method, gives in every world type of a fresh set of full families at least as many correct read outcomes as the
  current method and no more wrong ones, with reads calibrated on two real libraries".
- **FAIL** otherwise. No partial claim.

### Known before running (development; v5' on 14 spent runs, `fast/r13dev2/`)

**v5' passes 12 of 14:** sets 8-10 under both profiles, and rounds 10 and 11's six runs.
- **Round 12's two runs still fail**, by their O2 errors, which v5' does not address:
  - A119b: INTACT 3, CROSS 5, SILENT 2 against 0;
  - public: CROSS 2 against 0.
- **Against v5, the branch prefilter:**
  - **BRANCH in round 12:** correct +78 and +75, errors 136 -> 79 and 161 -> 99. DUX4's novel-exon copy is no longer lost in
    giant mixed clusters.
  - **Elsewhere:**
    - REMOVED moves by at most 16 reads (round 11 A119b: 1,845 -> 1,829, still far above current's 1,545);
    - BRANCH by up to 6;
    - 6 of the 14 runs are identical.
  - **O3's longest world:** at most 78 s, against 699 s before (round 12, public OUTSIDE).
- **Round 12's open O2 cases, with no fix in v5':**
  - copies DNA-identical over a read with a retained intron;
  - retained introns absorbed by a paralog's path as one gap-compressed edit;
  - margin-2 coincidences in an 80-copy array, which the genome's own alignment also misplaces.
- **Untested in development:** pseudogene and retrocopy families of this kind, and a 43-copy array.

### Reported beside the bars (not gating)

- fates per read class and error sources, **tabulated, every world, both profiles**;
- O1's drops (of how many silent copies), undecidable copies, new copies and twins;
- arm I's calls, true or false, every world;
- round 10's scoring and copy-level scoring;
- leave-one-family-out margins, **both profiles**;
- the cost: `timing.json`, the integrated steps apart from the current method's, and O3's peak memory.

### After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1311.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-17)

Both runs completed as pre-registered, with no amendment. Outputs are in `fast/heldout13/{A119b,public}/` (frozen read-only after
scoring), wall 417.8 s (A119b, which also built the local reference) and 325.8 s (public).

### Verdict: FAIL (A119b). Not claimed.

Bars from `e2e.json` (primary scoring):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 4,709 / 4,737 | 0 / **3** | 4,622 / 4,671 | 0 / 0 |
| CROSS | 4,783 / 4,860 | 2 / **3** | 4,690 / 4,781 | 3 / 0 |
| OUTSIDE | 3,990 / 4,383 | 683 / 16 | 3,817 / 4,230 | 786 / 8 |
| REMOVED | 4,657 / 4,682 | 10 / 0 | 4,467 / 4,618 | 95 / 1 |
| BRANCH | 4,718 / 4,888 | 95 / 63 | 4,622 / 4,785 | 102 / 72 |
| SILENT | 3,237 / 3,257 | 0 / **3** | 3,214 / 3,246 | 0 / 0 |

- **public passes all 6 worlds.**
- **A119b:**
  - **B1 holds in all 6 worlds,** by +20 to +393.
  - **B2 fails in three:** INTACT 3 vs 0, CROSS 3 vs 2, SILENT 3 vs 0.
- **Round 10's rule and copy level give the same verdict** (A119b FAIL, public PASS).

**The failure is one false O3 call, repeated in four worlds (found after the run, post hoc).** Every A119b error is arm I's PSV
call in MRGPRX on three reads of the present copy MRGPRX.1 (MRGPRX4). O2 makes no wrong-copy assignment in any of the 12 worlds.
- **The same three reads in every world.** The simulator draws an unchanged copy's reads from the same seed, so INTACT, CROSS,
  OUTSIDE and SILENT carry the same three reads under different names. No call is made in REMOVED, where MRGPRX.0 (whose window
  hosts the call) is the removed copy. BRANCH has the three reads too and makes no call; why was not traced.
- **They carry MRGPRX.1's heterozygous alleles.**
  - With a profile, `--het 0` means the profile's calibrated heterozygosity (`realworld.py`: `het_rate = 2 * het_per_base`),
    not a haploid donor.
  - MRGPRX.1's haplotype 2 carries two alleles 75 bp apart: C>G at chr11:18,269,598 and C>T at 18,269,673. 45 haplotype-2
    reads carry both, and 43 haplotype-1 reads carry neither.
- **Short reads turn two alleles into a delta-sized divergence.**
  - Of the 45 haplotype-2 reads covering both sites, O2 places 42 on MRGPRX.1. It leaves 3 unexplained: the two shortest (190
    and 231 bp) and a 506-bp read with 6 sequencing indels besides the two alleles.
  - Their primary alignments (MAPQ 60, on MRGPRX.1) have `de` 0.0158, 0.0129 and 0.0157. NM is 3, 3 and 11. That puts them in
    the world's top 0.2% (the median MAPQ-60 `de` is 0.0029, p99 0.0114).
  - **Both O2 stages leave them unexplained,** above delta plus floor: the paths and the genome fallback. Arm I clusters O2's
    unexplained reads, and without the call the three reads would be abstentions.
  - **The call:** arm I pools the family's unexplained reads across hosts; the three form a cluster on MRGPRX.0's window
    (`est_div` 0.024).
  - **Both checks pass:**
    - attribution keeps all three reads (its internals were not traced);
    - genome check: median `de` 0.0157 less the floor 0.0029 = 0.0128 > delta, and each read alone also passes the per-read check
      (the smallest is 0.0129 - 0.0029 = 0.0100).
  - The two alleles alone add 2/190 and 2/231 to `de` (0.0105 and 0.0087).
- **Delta is calibrated on whole copies.** It is the p99 exonic divergence between two haplotype assemblies over whole copies
  (register 1199). Two heterozygous sites within 200 bp of a short read exceed it, so the per-read genome check cannot separate an
  allele from a copy on short reads. This is the allele-versus-copy confound of the phased toy (2026-10-12). Round 10 met it too:
  TAF11L.33's haplotype-2 reads carried two donor alleles (median `de` 0.0103), and the error floor removed them (v4). Here the
  reads are short enough that two alleles exceed delta even after the floor.
- **Leave one family out:** A119b passes only without MRGPRX. Public passes with any family left out.
  - **Without OR7E (the 43-copy array), A119b's B1 margins turn negative:** INTACT -5, REMOVED -11, SILENT -11.
  - Its correct-read margin over the current method in those worlds is carried by one family.

**Reported beside the bars:**
- **Arm I:**
  - **REMOVED:** all 7 detectable removed copies called under both profiles (MRGPRX.0, YWHAZP.3, OR7E.15, SEPTIN7P.0, FTH1P.8,
    UBE2Q2P.3, CTSLP.1); PTMAP.8, within delta, not called.
  - **BRANCH:** all 6 novel-exon copies called under both profiles. PTMAP and FTH1P are single-exon retrocopies, with no BRANCH copy.
  - **Calls on present reads:** none under public. Under A119b, only the MRGPRX call above, in INTACT, CROSS, OUTSIDE and SILENT.
- **O1:**
  - **SILENT:** 24 of 27 silent copies dropped under both profiles; public keeps UBE2Q2P.1 as undecidable.
  - **OUTSIDE:** new copies in CTSLP, MRGPRX, OR7E, SEPTIN7P, UBE2Q2P and YWHAZP, with no twin.
  - **REMOVED:** new copies in FTH1P, SEPTIN7P (2), UBE2Q2P and YWHAZP.
- **Fates per read class** (correct / error / abstain / unknowable, current -> integrated):

| world | A119b catalogued | A119b other class | public catalogued | public other class |
|---|---|---|---|---|
| INTACT | 4709/0/499/0 -> 4737/3/468/0 | | 4622/0/479/0 -> 4671/0/430/0 | |
| CROSS | 4783/2/565/0 -> 4860/3/487/0 | | 4690/3/534/0 -> 4781/0/446/0 | |
| OUTSIDE | 3990/5/483/0 -> 4017/3/458/0 | uncatalogued 0/678/52/0 -> 366/13/351/0 | 3817/77/413/0 -> 3880/0/427/0 | uncatalogued 0/709/85/0 -> 350/8/436/0 |
| REMOVED | 3990/5/483/0 -> 4022/0/456/0 | absent detectable 636/5/43/0 -> 652/0/32/0; undetectable 31/0/2/13 -> 8/0/4/34 | 3817/77/413/0 -> 3886/0/421/0 | absent detectable 649/18/98/0 -> 727/1/37/0; undetectable 1/0/9/19 -> 5/0/4/20 |
| BRANCH | 4702/7/499/0 -> 4739/0/469/0 | absent detectable 16/88/119/0 -> 149/63/11/0 | 4622/0/479/0 -> 4671/0/430/0 | absent detectable 0/102/90/0 -> 114/72/6/0 |
| SILENT | 3237/0/311/0 -> 3257/3/288/0 | | 3214/0/304/0 -> 3246/0/272/0 | |

- **Cost: the branch prefilter holds on a 43-copy array.**
  - **The integrated steps, summed over six worlds:** 223 s (A119b) and 216 s (public), against round 12's 1,095 and 1,624 s.
    - O3: 164 and 160 s, longest world 72 s and 67 s (REMOVED), peak 0.17 GB.
    - O2 28-31 s (+ genome fallback 3-4 s); O1 22-27 s.
  - **The current method:** 188 and 183 s, peaking at 3.2 GB per world.
  - **Shared steps (simulate, map, mask):** 114 and 90 s. The local reference build (55 s, 14.4 GB) is paid once per set.

### Independent review (2026-10-17, from 01:04 PDT) and corrections

Review: `/mnt/linuxdisk/tmp/o3_graph/review13/REVIEW.md` (fresh subagent, read-only copies, its own scorer).

**The FAIL verdict stands.**
- Its scorer reproduces `e2e.json`, round 10's rule and copy level, both profiles: 0 differences over 305 cells.
- A119b fails in every scoring variant, in INTACT 3 vs 0, CROSS 3 vs 2 and SILENT 3 vs 0.
- The protocol is clean:
  - the pre-registration was committed while `heldout13` did not exist;
  - one run per profile, frozen code, no earlier exposure.
- Arm I re-runs byte-identical in all 6 A119b worlds.

**What the results said wrongly, also after the two corrections (a5d481a8, f4db7d5c):**
1. **"Two alleles exceed delta even after the floor" is wrong.**
   - The alleles alone are under every gate: 2/190 - 0.0029 = 0.0076 and 2/231 - 0.0029 = 0.0058.
   - Each read also needs 1, 1 or 6 sequencing indels.
   - Reverting the alleles on a scratch world makes two of the reads `assigned`.
2. **The call is a mosaic call, and its divergence is mostly one read's errors.**
   - The 3-read cluster has 1 unexplained alt (the PSV rule needs 2). The other allele recreates MRGPRX.0's base at a PSV,
     so the mosaic rule fires.
   - Its `est_div` (12/508 on MRGPRX.0's window, 18/506 on MRGPRX.2/.3's) counts the 2 alleles plus 10-16 positions covered
     only by the 506-bp read, where its sequencing indels count as the cluster's bases (`cluster_bases`, round 10's LCE1 rule).
   - Removing that read's indels removes the call.
3. **Attribution was traced.** The two alleles give exactly its fixed margin of 2 against MRGPRX.1. MRGPRX.2/.3's maps cover
   too few of the reads' positions to compete.
4. **The window and the pool are misdescribed.**
   - The call forms on MRGPRX.0's window and on MRGPRX.2/.3's. In OUTSIDE and SILENT it is reported on MRGPRX.2/.3's.
   - On MRGPRX.1's own window, the 3 reads sit in a 10-read mixed cluster.
   - It depends on the family's pool:
     - dropping any one of three unrelated MRGPRX.0 reads removes it;
     - BRANCH minus its 6 MRGPRX.B reads recreates it;
     - REMOVED with INTACT's pool composition recreates it on MRGPRX.2/.3's window.
   - The REMOVED explanation in the results is refuted.
5. **"No O2 error in any world" and "every A119b error is arm I's" hold only for catalogued reads and the failing cells.** Integrated
   O2 places OUTSIDE's uncatalogued within-delta PTMAP.8 reads on PTMAP.7 (13 / 8). In BRANCH it places novel-exon-free reads
   of the branch copies on their parents (63 / 72). Both are DNA-identical placements.
6. **Round 10's comparison is incomplete.** Round 10's TAF11L.33 call was also a present copy's haplotype-2 reads called on
   another copy's window; its fix was attribution's no-rival rule (88c445d9). That shape is the recurring mechanism.

**Also found:**
- **Public's PASS cannot test the allele confound.** `het_sites` draws `round(2 x het_per_base x exonic bp)` sites per copy.
  Public's rate rounds to 0 for 86 of 89 copies, and none has two. A119b gives 7 copies with two or more.
- **Arm I pools O2's paths-stage table.**
  - Its path aligner (`minimap2 -x map-hifi`, minimum score 200) cannot align reads under ~200 bp.
  - 385 of the 436 path-unexplained reads in A119b INTACT have no path alignment.
  - The genome fallback leaves 42 unexplained. Arm I's pool is therefore about 10x what O2 finally leaves.
- **Per-family losses under the B1 margins:**
  - A119b INTACT: FTH1P -7, MRGPRX -3, PTMAP -3;
  - A119b REMOVED: CTSLP -16. CTSLP.1's undetectable reads are placed by the integrated pipeline (`unknowable`), while the
    current method's call holds them.
- **Zero OUTSIDE errors on FTH1P.8 rest on the genome check.** Arm I's unchecked calls hold 312 (A119b) and 396 (public)
  uncatalogued reads.

**Fix directions** (to be developed on all spent runs, round 13 included, and pre-registered on a fresh set):
- the mosaic distance should not take one read's bases as the cluster's;
- attribution's margin should grow with read length (the existing `--member-attribution-delta`);
- pool composition: pool only what the genome fallback leaves unexplained, or require a call to hold on its members' own host
  window;
- an allele-aware check is right in principle but not as the only lever: it can veto true calls (v7 lost CTSLP.1), and it
  depends on O2's placements;
- fix the simulator's het-site rounding, pre-registered, so that the public donor is not effectively haploid.


---

## PREREG_o3_graph_round14_2026-10-17.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round14_2026-10-17.md`)*

## Round 14: the integrated pipeline v9 on a fifteenth set of full families, with a diploid donor (2026-10-17)

Pre-registration, written before any world of the fifteenth set was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** Round 13 (register 1311) failed under A119b on one false O3 call. Its review traced the call to O3's mosaic distance,
which counted one long read's sequencing indels at positions only that read covers, on a present copy's haplotype-2 reads.
Round 12 (register 1310) failed on O2 errors: retained-intron reads that a paralog's path absorbed as one gap-compressed edit,
and two coincidences in an 80-copy array. v9 adds one rule against each mechanism that development could fix. The simulator's
donor is now diploid under both profiles, and development under that donor found the allele confound again, which v9's
haplotype check addresses.

**Question.** On a fresh set of full families, does v9 do at least as well as the current method in every world type, with both
calibrated v2 profiles: at least as many correct read outcomes and no more wrong ones?

### The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v9** = round 13's v5' + three rules:
- **O3 `--mosaic-shared`** (551c7685). The mosaic rule's nearest-copy distance is taken only over positions that >= 2 of the
  cluster's members cover. One read's sequencing errors where it alone covers are no evidence of the cluster's copy.
- **O3 `--haplotype-check`** (1c604e47 + cd93a327). A PSV call is dropped when the haplotype it is stays within delta of a copy:
  - the family's reads that O2 explains are aligned to the call's host window;
  - those carrying one of the call's unexplained alts (alts no copy's map carries) join it;
  - the extended cluster's nearest-copy distance, over positions >= 2 of them cover, must exceed delta.
  A present copy's haplotype-2 reads carry its donor alleles over the whole copy, whereas a missing copy's PSVs are carried
  by its own reads.
- **O2 `--o2-genome veto`** (fe6804cc; `o2graph.py run W integv 2`, replacing the fallback step `integl`). It is the fallback
  for reads the paths leave unexplained, and it also rechecks every read the paths assign, on the read's own genome
  alignments. A path call is overruled only when three things hold:
  - another seeded copy's locus, or a locus off them, beats the called copy's by the margin, within the gate; clipped ends and
    long insertions cost their length (`veto_distance`);
  - the read's primary alignment lies on that winner;
  - at MAPQ >= 30.
  A genome tie or a thin margin keeps the path call.
- Everything else as v5':
  - `o1seed.py --twins`;
  - `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`;
  - `detect.py --member-genome-check --member-attribution --error-floor --o1-hosts --branch-genome-novel --arms I`.
- **Measured and not adopted** (each turned a passing development run into a failure, set 10 A119b REMOVED):
  - `--allele-discount` v7f;
  - `--member-attribution-delta`;
  - `--o3-pool final` (it does not remove round 13's call, and it adds REMOVED and BRANCH errors).

**The simulator (both systems): `--sim-het-draw`** (551c7685). Each copy's het and hom site counts are drawn (Binomial over its
exonic bases) instead of rounded. Round 13's review found the public donor effectively haploid under rounding: 86 of 89 copies
had no het site. The calibrated rates are unchanged.

```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared
    --haplotype-check --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```

### Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout14`, e69d201e), all in no earlier
list:
- segmental-duplication pseudogene families: SMG1P, PKD1P, CHEK2P, RRN3P, CDC27P;
- few-copy gene families: AKR1C, HLA-DQ, PCDHB, ANKRD36, TPTE, NXF, SERPINA, MICAB;
- olfactory subfamilies: OR51, OR8, OR10, OR5;
- retrocopy families as spares: NPM1P, HMGB1P, GAPDHP, RPL21P, HNRNPA1P.

**Build rule:** round 12's, `families.py build families14_all.json --cap 1000 --link locus --min-cover 0.8`, every group (102 s).
The set is the first 8 families with >= 3 copies, in list order:

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| SMG1P | 5 | 14 | chr16 |
| PKD1P | 3 | 6 | chr16 |
| CHEK2P | 7 | 11 | chr2, chr10, chr15, chr22, chrY |
| CDC27P | 13 | 13 | chr2, chr13, chr14, chr15, chr20, chr21, chr22, chrY |
| PCDHB | 14 | 19 | chr5 |
| OR8 | 3 | 5 | chr11 |
| OR10 | 4 | 5 | chr11 |
| NPM1P | 25 | 25 | 13 chromosomes |

74 copies. RRN3P, AKR1C and OR5 gave 2 copies; HLA-DQ, ANKRD36, TPTE, NXF, SERPINA, MICAB and OR51 no family.
Files: `/mnt/linuxdisk/tmp/o3_graph/families14.json` (this set), `families14_all.json`, `families14.log`.

### Worlds

Six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
`profiles/public_testis_v2.json` (`--het 0`: the profile's calibrated heterozygosity): 12 worlds.

**The seed:** `precheck.py families14.json heldout14 201 211` (`precheck14.txt`; 96 s, 13.3 GB) gives the smallest seed in
201-210 with >= 3 detectable removed copies: **seed 201**.
- Detectable: OR8 0.1158, OR10 0.0834, NPM1P 0.0529, PCDHB 0.0324, CHEK2P 0.0161.
- Within delta: PKD1P 0.0019, SMG1P 0.0003, CDC27P 0.0.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families14.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout14
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout14 201 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout14 201 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```

### Read outcomes, bars and reading

Round 12's, unchanged. The e2e scorer reads the veto's table (`--o2 o2.integv.m2.tsv`).
- **Read outcomes:** per-read detectability less the world's error floor. Round 10's rule and copy level are secondary.
- **Bars:** B1 (integrated correct >= current) and B2 (integrated errors <= current), per world, all 12 worlds.
- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "the integrated O1 -> O2 -> O3 pipeline v9, which uses nothing from the
  current method, gives in every world type of a fresh set of full families at least as many correct read outcomes as the
  current method and no more wrong ones, with reads calibrated on two real libraries".
- **FAIL** otherwise. No partial claim.

### Known before running (development)

**Development data:** 16 spent runs (sets 8-10 and rounds 10-13, both profiles, two untrimmed-library runs; round 13 is
development since its review), plus 8 spent family sets re-simulated with the diploid donor (`--sim-het-draw`: rounds 9-11
and 13, both profiles). Each lever was measured alone against v5'. The adoption rule was fixed before any result was read: a
lever must remove round 13's false call or reduce errors, and must not turn a passing run into a failure
(`docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 14, `fast/r14dev/`).

| | 16 spent runs passing | het-draw runs passing |
|---|---|---|
| v5' (round 13) | 13 | not run |
| v8 = v5' + `--mosaic-shared` + veto | 15 | round 13's set: A119b FAIL (FTH1P.6 haplotype), public PASS |
| **v9 = v8 + `--haplotype-check`** | **15** | **8 of 8** |

- **The one development failure is round 12 A119b.** Its 2 DUX4 margin-2 coincidences and 1 KRTAP2 read, DNA-identical at
  two copies, are placed by O2 in INTACT, CROSS and SILENT. The current method abstains on them, and no sequence rule
  separates them.
- **v9's O3 rules against v5' on the 16 runs** (O2 unchanged): 13 runs are identical and 3 change, all for the better:
  - round 13 A119b: its false call removed, REMOVED +2;
  - round 11 A119b: REMOVED +11;
  - set 10 A119b: REMOVED +45.
- **The veto:**
  - it fixes round 12's four CDRT15 retained-intron errors (round 12 public passes);
  - it cuts OUTSIDE errors in 8 runs;
  - it costs 1 and 2 CROSS reads in set 8.
- **The haplotype check took three versions:**
  - in genome coordinates it missed a 45-copy array (TAF11L.10's haplotype-2 reads, whose MAPQ-0 primaries lie on three
    units);
  - selecting carriers by every alt pulled in other copies' reads;
  - the adopted version selects by unexplained alts on the host window.
- **Not adopted:**
  - `--allele-discount` and `--member-attribution-delta` each turned set 10 A119b into a failure;
  - `--o3-pool final` kept round 13's call and added errors.
- **Cost:** the haplotype check aligns the family's explained reads to each host window that has a candidate call. On round
  13 A119b, O3's longest world went from 26 to 58 s.
- **Untested in development:** these families; segmental-duplication pseudogene families; a protocadherin cluster; the diploid
  donor on round 12's and sets 8's families.

### Reported beside the bars (not gating)

- fates per read class and error sources, tabulated, every world, both profiles;
- O1's drops (of how many silent copies), undecidable copies, new copies and twins;
- arm I's calls, true or false, every world;
- the veto's overrulings and the haplotype check's removals, right or wrong;
- round 10's scoring and copy-level scoring;
- leave-one-family-out margins, both profiles;
- the cost: `timing.json`, the integrated steps apart from the current method's, and O3's peak memory.

### After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1312.
- Results appended below a line, without editing anything above it.

---

## Results (2026-10-17)

Both runs completed as pre-registered, with no amendment. Outputs are in `fast/heldout14/{A119b,public}/`, frozen read-only
after scoring. Wall time was 305 s for A119b, which also built the local reference, and 232 s for public.

### Verdict: PASS under both profiles (12 of 12 worlds). Claimed.

Bars from `e2e.json` (primary scoring):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 3,389 / 3,447 | 0 / 0 | 3,334 / 3,391 | 0 / 0 |
| CROSS | 3,449 / 3,513 | 0 / 0 | 3,387 / 3,451 | 0 / 0 |
| OUTSIDE | 3,308 / 3,588 | 232 / 39 | 3,288 / 3,574 | 220 / 32 |
| REMOVED | 3,502 / 3,584 | 0 / 0 | 3,472 / 3,558 | 4 / 0 |
| BRANCH | 3,393 / 3,691 | 248 / 125 | 3,339 / 3,610 | 270 / 143 |
| SILENT | 2,551 / 2,609 | 0 / 0 | 2,584 / 2,636 | 0 / 0 |

- **B1 holds in all 12 worlds,** by +52 to +298.
- **B2 holds in all 12.** In INTACT, CROSS, SILENT and REMOVED the integrated pipeline makes no error under either profile.
- **The donor is diploid under both profiles.** 51 of A119b's 74 expressed copies and 19 of public's have haplotype-2 reads
  (round 13 under rounding: public 3 of 89).

**Reported beside the bars:**
- **Secondary scorings:**
  - public passes under round 10's rule and at copy level;
  - **A119b fails both in REMOVED B2:** 1 vs 0 errors (round 10's rule) and 3 vs 0 (copy level). Every other cell holds.
- **Leave one family out: both profiles pass in all 16 rows.**
  - CHEK2P carries most of A119b's margins. Without it: REMOVED +0 (B1 holds with equality), INTACT +3, CROSS +6, SILENT +2.
  - Without it, public is at INTACT +6 and SILENT +2.
- **Where the integrated errors are:**
  - OUTSIDE: 39 and 32 reads of uncatalogued copies placed on a catalogued paralog;
  - BRANCH: 123-142 reads of branch copies that do not carry the novel exon, placed on the parent copy, plus 1-2 placed by
    the genome step.
- **Arm I:**
  - **REMOVED:** every detectable removed copy called under both profiles (CHEK2P.4, PCDHB.0, OR8.2, OR10.0, NPM1P.12). No
    call holds a present read in any world.
  - **BRANCH:** all 7 novel-exon copies called (SMG1P, PKD1P, CHEK2P, PCDHB, OR8, OR10, NPM1P). CDC27P has no branch copy.
- **O1:**
  - **SILENT:** 18 of 22 silent copies dropped; CHEK2P.2, PCDHB.1 and PKD1P.0 kept as undecidable.
  - **OUTSIDE:** new copies in 6 families (A119b) and 7 (public), one of them a twin (SMG1P).
- **Fates per read class** (correct / error / abstain / unknowable, current -> integrated):

| world | A119b catalogued | A119b other class | public catalogued | public other class |
|---|---|---|---|---|
| INTACT | 3389/0/514/0 -> 3447/0/456/0 | | 3334/0/505/0 -> 3391/0/448/0 | |
| CROSS | 3449/0/551/0 -> 3513/0/487/0 | | 3387/0/536/0 -> 3451/0/472/0 | |
| OUTSIDE | 3308/0/282/0 -> 3350/0/240/0 | uncatalogued 0/232/81/0 -> 238/39/36/0 | 3288/1/264/0 -> 3334/0/219/0 | uncatalogued 0/219/67/0 -> 240/32/14/0 |
| REMOVED | 3314/0/276/0 -> 3367/0/223/0 | absent detectable 180/0/48/0 -> 217/0/11/0; undetectable 8/0/33/44 -> 0/0/16/69 | 3294/0/259/0 -> 3342/0/211/0 | absent detectable 178/4/38/0 -> 216/0/4/0; undetectable 0/0/24/42 -> 0/0/0/66 |
| BRANCH | 3389/0/514/0 -> 3447/0/456/0 | absent detectable 4/248/139/0 -> 244/125/22/0 | 3334/0/505/0 -> 3391/0/448/0 | absent detectable 5/270/103/0 -> 219/143/16/0 |
| SILENT | 2551/0/464/0 -> 2609/0/406/0 | | 2584/0/467/0 -> 2636/0/415/0 | |

- **Cost:**
  - **The integrated steps, summed over six worlds:** 159 s (A119b) and 149 s (public).
    - O3: 98 and 100 s, longest world 44 s, peak 0.28 GB.
    - O2 18-21 s (+ the genome veto 18-23 s); O1 11-19 s.
  - **The current method:** 133 and 129 s, peaking at 2.5 GB per world.

### Independent review (2026-10-17, from 09:30 PDT) and qualifications

Review: `/mnt/linuxdisk/tmp/o3_graph/review14/REVIEW.md` (fresh subagent, read-only copies, its own scorer, end-to-end
reproductions).

**The PASS stands as a pre-registered result, but it is one fragile draw.**
- **What holds:**
  - the reviewer's scorer reproduces all three scorings, both profiles, cell for cell (0 differences in 268 cells);
  - the protocol is clean: frozen code, one launch with the pre-registered flags, no re-runs, no trace of these families in
    any earlier truth file;
  - both runs reproduce byte-identically end to end, also with the current method's outputs hidden.
- **Post-hoc seed probe** (frozen code, same families, seeds 202-205; it cannot change the verdict): public passes 4 of 4.
  **A119b passes 2 of 4:**
  - seed 202 fails REMOVED 6 vs 0: present NPM1P.14 reads with drawn homozygous alleles, called on the window of a new copy
    that O1 seeds at the NPM1 parent gene; the genome check passes them by +0.0005;
  - seed 205 fails INTACT and CROSS 4 vs 0: CDC27P.8's haplotype-2 reads called on a paralog's window; the haplotype check
    misses by 0.00024 (0.00982 vs delta).
  Both are round 13's mechanism, a present copy's donor-allele reads called on another locus's window. **The allele confound is
  not solved.**
- **What decided this round's cells:**
  - public CROSS rests on `--haplotype-check` alone: without it B2 fails 4 vs 0, so round 13's v5' would have failed. The
    removed call is 4 hap-1 reads of CHEK2P.5, a 246-bp copy, carrying a homozygous allele that `--sim-het-draw` drew; the
    genome check let it through by +0.0011;
  - the veto overruled 0 reads in 12 worlds (`o2.integv` = `o2.integl`), and `--mosaic-shared` changed no bar: round 14 is no
    evidence for either;
  - A119b REMOVED B2 (0 vs 0) is decided by one placed PCDHB.0 read 0.00008 under the detectability threshold. At floor 0.0031
    instead of 0.0032 the cell fails. Round 10's rule and copy level fail it (1 and 3 vs 0);
  - most of the B1 margin in INTACT, CROSS, SILENT and REMOVED is one copy's reads (CHEK2P.5, which copy_assign mostly calls
    ambiguous). Without CHEK2P the margins are +0 to +14.
- **Corrections to the results:**
  - BRANCH's integrated errors are not DNA-identical: 105 of 125 (A119b) and 127 of 143 (public) carry the branch copy's SNPs.
    OUTSIDE's are DNA-identical.
  - The pre-registered report of the veto's overrulings and the haplotype check's removals is missing. It is given above: 0
    overrulings; removals = the public CROSS call above, plus PCDHB.B PSV calls whose reads are also in PCDHB's branch call.
- **Undisclosed in the pre-registration:** the family set and the seed were fixed (01:33-01:40) 7.5 h before the method was
  frozen, during v8/v9 development. The reviewer found no use of them in any development run.

**The claim is literally true for this draw under the primary scoring and must be cited with these qualifiers:**
- one simulated draw; at 4 further seeds, A119b fails 2;
- the two secondary scorings fail A119b REMOVED;
- public CROSS depends on the haplotype check;
- the veto and the shared-position mosaic distance were not exercised;
- it is no evidence that the allele confound is solved.


---

## PREREG_o3_graph_round15_2026-10-17.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round15_2026-10-17.md`)*

## Round 15: integrated v11 on fresh families, three seeds per profile, and the cost of full families (2026-10-17)

Pre-registration, written before any world of round 15's families was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.**
- **Round 14 passed (register 1312), but its review called the pass one fragile draw:** A119b passed only 2 of 4 further
  seeds. v10 then passed all 18 spent seed and diploid-donor runs.
- **The goal check found the cost claim unshown on full families.** On the NPIP-sized mainstay set, the integrated stages took
  1,074 s against the current method's 369 s (A119b).
- **Round 15 tests both claims, each with its own pre-registered bar.**

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v11** = v10 + two changes. v10 is round 14's v9 plus O3 `--haplotype-check-genome` (f867ec7c). Both changes were
developed on spent runs only:
- **O3 `--branch-linear`** (44798db2): the branch step without the family graph.
  - Candidates are spliced-aligned (minimap2 splice:hq) to the family's loci as linear sequences.
  - A read base is novel only where no locus alignment explains it, which stands for the graph's path switches between copies.
  - vg msga is gone. It was 121 of O3's 158 s on the mainstay CROSS world, and it was not deterministic.
  - **On all 34 spent runs** (16 rounding runs, round 14's 10 seed runs, 8 diploid-donor sets): every verdict is v10's, 14
    runs have identical bars, and 19 gain 1-9 correct reads with equal or fewer errors.
- **`fastset --aln-cache`** (f2e59107, 7839138c): exact alignment caches shared by a set's worlds.
  - O2's path alignments are cached per read sequence.
  - O3's catalog-isoform-to-host alignments are cached by window, query names and sequences.
  - The O2 tables were byte-identical, cold and warm, on the mainstay worlds.

**Measured on spent runs and not adopted**, because each turned a passing development run into a failure:
- **O2 `--twin-abstain`** (de03bf5e): a read that neither of two twins beats by the margin is abstained. Round 10 public
  OUTSIDE fell to 751 correct, below the current method.
- **O2 `--o2-genome guard`** (de03bf5e): a kept path call is abstained when a locus off every seeded copy explains the read as
  well. Round 10 A119b and public fell below the current method's correct count in several worlds; the guard abstains on
  silent array units' ties.

Both cut the O2 proof's OUTSIDE errors on near-identical twins, from 138-171 to 0-2 with the two together. That failure stays
open: an identical, uncatalogued, expressed twin gives O1 no uniquely placed read to seed it.

**Integrated flags:**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

### Families (held out; drawn by a rule fixed before the pool was built)

**The pool** (`candidates15.py`, e1fd5bb5):
- Every CAT protein-coding gene-name stem with >= 3 genes, none of them matched by any earlier group: 1,409 stems. The list's
  sha1 is 808e8f6d.
- Built with round 12's rule, `families.py build --cap 1000 --link locus --min-cover 0.8` (543 s).
- 105 stems formed a family; 30 have >= 3 copies (eligible), 7 of them >= 6 copies.

**The draw:**
- **Rule** (fixed before the build, applied once): eligible names sorted; seed 15. When >= 4 eligible families have >= 6
  copies, draw 4 of those, then 4 of the rest.
- **Large families drawn:** IFNA, FAM90A, SCYGR, SSU72P, from FAM90A, IFNA, OR2L, OR4A, SCYGR, SSU72P and VN1R.
- **Others drawn:** H2AB, KRT, ANP, BOLA, from 23.

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| IFNA | 7 | 7 | chr9 |
| FAM90A | 57 | 65 | chr8, chr12 |
| SCYGR | 10 | 10 | chr2 |
| SSU72P | 7 | 7 | chr7, chr11 |
| H2AB | 3 | 3 | chrX |
| KRT | 4 | 4 | chr17 |
| ANP | 3 | 6 | chr4, chr12, chr15 |
| BOLA | 3 | 12 | chr16 |

- 94 copies. FAM90A is the family the user parked for O2/O3 on 09-17, as a good near-identical copy-assignment case.
- **Files:** `/mnt/linuxdisk/tmp/o3_graph/families15.json` (this set) and `r15/families15_all.json`, `r15/build.log`.

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`.
- **Seeds:** 501, 502 and 503, fixed with no seed selection. That makes 6 cells (seed x profile) and 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families15.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout15
for S in 501 502 503; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r15 $S --local $O/local $V --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

### Bars

**Correctness.** Read outcomes are round 14's: per-read detectability less each world's error floor, scored on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 15 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1).** Measured on the mainstay set: NPIP as one family, TBC1D3, GOLGA8, GOLGA6L, CT47A, XAGE and SPATA31; seed
306; both profiles. These families are spent, but cost is not tuned on them.
- **Setup:** `fastset.py run` with v11's flags, the existing mainstay local reference, nothing else running.
- **Bar:** v11's integrated wall time, summed over the six worlds (O1 + O2 + O2_genome + O3), is at most the current method's
  summed wall time, in each profile.
- **Reported, no bar:** v10's run on the same set, and every step's wall time and peak RSS (`timing.json`).

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that v11 does at least as well as the current method in every world type, on 8 held-out
families, under both profiles and three seeds.

**A cost PASS** would claim that, on NPIP-sized families, the integrated O1 -> O2 -> O3 costs no more wall time than the
current method on the same machine and settings.

**Neither would claim:**
- anything on real reads;
- that near-identical uncatalogued twins are handled (see above).

### Amendment 1 (2026-10-17 17:25, before any rerun)

**What happened:** /mnt/linuxdisk filled up during the batch (0 bytes free). Cell A119b 501 finished; every later run, cells
2-6 and the v10 cost runs, exited with code 120 within a second, at its first write. C1's two v11 runs had already finished.

**The space freed:** 37 per-run masked genome copies (`masked/masked.fa`, 3.1 GB each, 108 GB in all) of runs completed and
recorded on 10-15 and 10-16, none referenced by an active directory.
- Each copy is the CHM13 reference with a few intervals hard-masked. It can be regenerated exactly by `realworld.py mask` from
  the run's families, set and seed.
- The list of deleted copies is `docs/deleted_masked_2026-10-17.txt`.

**Under the rule above:** the crashed runs are rerun once, unchanged. Cell public 501's partial directory is removed first:
its mask step died writing.

### Results (2026-10-17; cells 2-6 rerun once under Amendment 1)

**Verdict: FAIL.** 3 of the 6 cells pass; the bar was all six. **C1 FAILS too.**

| cell | verdict | failing worlds (integrated errors vs current) | integrated / current wall |
|---|---|---|---|
| A119b 501 | PASS | | 855 / 417 s |
| public 501 | PASS | | 662 / 374 s |
| A119b 502 | FAIL | CROSS 3 vs 0, INTACT 2 vs 0 | 944 / 468 s |
| public 502 | PASS | | 771 / 411 s |
| A119b 503 | FAIL | CROSS 5 vs 0, INTACT 5 vs 0, SILENT 5 vs 0 | 1,117 / 462 s |
| public 503 | FAIL | CROSS 2 vs 0, INTACT 1 vs 0, SILENT 1 vs 0 | 934 / 411 s |

**B1 holds in all 36 worlds**: integrated correct exceeds current everywhere.

**Every failing bar is a few errors where the current method makes none, and every one of those errors is in FAM90A** (57
near-identical copies). By the O2 proof's oracle:
- **Most are undecidable reads.** The read's alignment cost ties at the true copy and a sibling, yet path O2 assigned it by the
  margin, to a copy with a richer path set:
  - A119b 503: five reads of copies 20, 21, 29 and 41, all assigned to FAM90A.45;
  - public 503: copy 29 -> 39, copy 12 -> 11.
- **One read is decidable** but was called against the oracle (A119b 502: FAM90A.45 -> FAM90A.2).
- **This is a calibration gap in large near-identical families.** Path structure asymmetries that union projection does not
  remove decide the margin, and the margin is fixed (2 edits) whatever the number of rivals. copy_assign's PSV certificate
  divides alpha by the number of rivals, and it abstains on these reads.
- **Elsewhere the integrated pipeline stays far ahead:** OUTSIDE 4-66 errors against 101-263; BRANCH 13-60 against 44-195.

**C1, the cost on the mainstay set:** v11's integrated stages took 416 s (A119b) and 418 s (public), against the current
method's 320 s and 333 s.
- **That fails C1** by 96 s and 85 s.
- **It is 2.4-2.9x faster than v9's integrated stages** (1,200 s and 987 s): O3 825 -> 149 s and 763 -> 98 s, without the vg
  graph.
- **O2 is now the larger cost** (228 / 279 s). The cache helps only worlds whose path set an earlier world made.
- **The thread counts were not equal:** the current method's minimap2 ran 8 threads under `--heavy-current`, and each integrated
  step had 2. Both runs' bars are identical to v9's (PASS).
- **On round 15's families,** where FAM90A's 57 copies make O2 large, the integrated stages took 1.7-2.4x the current method.
- **Reported, no bar: v10 on the mainstay set** (rerun under Amendment 1, verdict PASS in both profiles). Integrated stages
  1,057 s (A119b) and 1,026 s (public), current 306 s and 313 s; O3 alone 816 s and 746 s, which `--branch-linear` cut to
  149 s and 98 s in v11.


---

## PREREG_o3_graph_round16_2026-10-17.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round16_2026-10-17.md`)*

## Round 16: integrated v12 (the O2 margin) on every unused family of round 15's pool, and the cost of full families in CPU seconds (2026-10-17)

Pre-registration, written before any world of round 16's families was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.**
- **Round 15 failed 3 of 6 cells (register 1316).** Every failing bar was a few FAM90A errors where the current method makes
  none. Read by read, path O2's fixed margin of 2 edits:
  1. rewarded a copy for read-derived paths its siblings lacked;
  2. took a locus-fallback call on bases its own alignment clipped;
  3. took a call that no rival's paths could answer.
- **Round 12 A119b, the one development run v10 and v11 never passed,** adds two more mechanisms:
  4. a margin of two 1-bp homopolymer indels against 63 DUX4 rivals;
  5. a retained intron that only one copy's annotation absorbs.
- **Four rules in O2's genome step address them.** They were developed on spent runs only (section 16 of
  `docs/O1O2O3_FAST_REALISTIC_2026-10-16.md`).
- **Round 15's cost bar compared wall times at unequal threads:** the current method's minimap2 ran at 8 threads, each
  integrated step at 2. Round 16's cost bar is in CPU seconds, which `fastset` now records for every step (d6a8e6ba).

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v12** = v11 + four O2 rules. Each acts in the genome step (`o2.integv`), so O3's input (`o2.integ`) is unchanged:
- **`--veto-abstain`** (1bd6693a, mechanism 1): a path call is abstained when a rival's worst-case genome cost (every clipped
  base paid) beats the call's best case (clips under 30 bases free) by the margin, and no confident winner takes the read.
- **`--locus-pairwise`** (10208f8c, mechanism 2): against each rival, a locus-fallback call pays its edits plus the read bases
  the rival aligns and it clips; it must still win by the margin.
- **`--o2-significance`** (1bd6693a, mechanism 4): each path edit is priced at -ln of its type's rate in the world's uniquely
  placed reads (mismatch, homopolymer indel, other indel). A call must beat every rival by ln(rivals / alpha), alpha = 1e-3,
  as copy_assign's certificate divides alpha by its rivals.
- **`--structure-check`** (10208f8c, mechanisms 3 and 5): a rival whose own locus explains the read better than its paths
  do, by the margin, is compared with the call on the genome; the call is abstained unless it wins there by the margin.

**On all 46 spent development runs (276 worlds), against v11's table** (register 1317):
- 36 -> 40 runs pass. The four that turned are round 15 A119b 502 and 503, round 15 public 503, and round 12 A119b. No world
  that passed fails.
- Correct -1,294 (-0.21%) and errors -139 (-2.1%), of v11's 603,857 and 6,617.
- The six still failing are the O2 proof's cells, on their OUTSIDE twins.
- The thinnest B1 margins left are 3-6 reads.

**Integrated flags (`$V12` below):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check
--sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

### Families (held out; no draw)

**Every family of round 15's pool with >= 3 copies that round 15 did not draw** (`pool16.py`, d6a8e6ba):
- **The pool** is round 15's (`r15/families15_all.json`): every unused CAT protein-coding stem with >= 3 genes, built by
  round 12's rule. 105 stems formed a family and 30 have >= 3 copies; round 15 drew 8 of those.
- **Round 16 takes the other 22.** Taking all of them leaves nothing to select.
- **22 families, 90 copies.** The largest are OR4A (7 copies), OR2L and VN1R (6). FAM90A- and DUX4-scale arrays are not in
  this set; they remain tested on spent runs only.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families16.json` (set `r16`, sha1 3698cae80929).

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| CLEC | 5 | 20 | chr16 |
| CST | 3 | 4 | chr20 |
| GGTLC | 4 | 8 | chr20, chr22 |
| KCNJ | 3 | 4 | chr17 |
| LGALS | 4 | 16 | chr17 |
| MBD3L | 5 | 5 | chr19 |
| MKRN | 3 | 6 | chr7, chr12, chrX |
| MTRNR2L | 4 | 4 | chr3, chr5, chr11, chr17 |
| MYL | 3 | 9 | chr6, chr18 |
| NUTM | 4 | 10 | chr10 |
| OR13C | 3 | 3 | chr9 |
| OR2J | 4 | 5 | chr6 |
| OR2L | 6 | 9 | chr1 |
| OR2M | 5 | 6 | chr1 |
| OR4A | 7 | 7 | chr11 |
| OR4N | 3 | 3 | chr14, chr15 |
| SLC35G | 4 | 4 | chr8, chr17, chr18 |
| SLC6A | 3 | 10 | chr16, chrX |
| SPAG | 3 | 12 | chr8 |
| TAS2R | 4 | 4 | chr12 |
| TRPC | 4 | 4 | chrX, chrY |
| VN1R | 6 | 6 | chr13, chr21, chr22 |

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`.
- **Seeds:** 601, 602 and 603, fixed with no seed selection: 6 cells (seed x profile), 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families16.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout16
for S in 601 602 603; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r16 $S --local $O/local $V12 --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

### Bars

**Correctness (as round 15).** Read outcomes are per-read detectability less each world's error floor, scored on the veto
table (`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 16 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds.** Measured on the mainstay set (NPIP as one family, TBC1D3, GOLGA8, GOLGA6L, CT47A, XAGE,
SPATA31; seed 306; both profiles), as in round 15. These families are spent, but cost is not tuned on them.
- **Setup:** `fastset.py run` with v12's flags, the existing mainstay local reference, nothing else running:
  ```
  M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
  for P in A119b public; do
    python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r16cost/v12/$P \
      /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V12 \
      --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
  done
  ```
- **The measure:** each step's `cpu` in `timing.json` (fastset's `/usr/bin/time` %U + %S).
- **Bar:** v12's integrated CPU seconds (user + system, children included), summed over the six worlds (O1 + O2 +
  O2_genome + O3), at most the current method's summed CPU seconds, in each profile.
- **Reported, no bar:** wall times (at the flags' threads: current 8, each integrated step 2), peak RSS, and the CPU seconds
  of round 16's six cells.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that v12 does at least as well as the current method in every world type, on 22
held-out families, under both profiles and three seeds.

**A cost PASS** would claim that, on NPIP-sized families, the integrated O1 -> O2 -> O3 uses no more CPU than the current
method.

**Neither would claim:**
- anything on real reads;
- that near-identical uncatalogued twins are handled (the O2 proof's end-to-end failure, register 1314-1315);
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs.

### Results (2026-10-17, 20:33-21:16; nothing changed after the first run, no crash)

**Verdict: FAIL.** 4 of the 6 cells pass; the bar was all six. **C1 PASSES** in both profiles.

| cell | verdict | failing worlds (integrated vs current correct) | integrated / current CPU |
|---|---|---|---|
| A119b 601 | FAIL | INTACT 4,232 vs 4,238; SILENT 3,224 vs 3,231 | 148 / 240 s |
| public 601 | PASS | | 119 / 211 s |
| A119b 602 | FAIL | SILENT 2,638 vs 2,642 | 168 / 287 s |
| public 602 | PASS | | 132 / 270 s |
| A119b 603 | PASS | | 167 / 273 s |
| public 603 | PASS | | 135 / 245 s |

**B2 holds in all 36 worlds.** The integrated pipeline makes no error where the current method makes none, and far fewer
elsewhere: OUTSIDE 44-182 against 811-1,214, REMOVED 0-3 against 39-132, BRANCH 157-304 against 285-653. Over the 36 worlds it
makes 1,833 errors against 8,796, and 147,583 correct calls against 137,426 (+7.4%).

**Every failing bar is B1, by 4-7 correct reads, in INTACT or SILENT under A119b.** Post-hoc diagnosis (read-only, on the
round's tables):
- **The two methods nearly tie there.** In A119b 601 INTACT, 28 reads are correct under the current method only and 22 under
  the integrated pipeline only.
- **Most of the 28 are reads integrated O2 leaves `unexplained`:** more edits than the fixed divergence gate,
  (delta + floor) x length, allows. The rest are `thin` (margin 1) or locus ties. Round 10's set 10 failure was the same gate.
- **v12's new rules cost 5 of the 13 missing reads (`structure_only`), but not the verdict.** Rebuilt with v11's genome
  step, the cell still fails: INTACT 4,235, SILENT 3,226.
- **The integrated lead is thin in these worlds because the families are small** (3-7 copies). The current method's
  certificate assigns single-PSV reads among few rivals, and short noisy reads.

**C1, the cost on the mainstay set** (nothing else running): v12's integrated stages used 810 CPU seconds (A119b) and 660
(public), against the current method's 1,046 and 949.
- **That passes C1:** 23% and 30% less CPU.
- **Wall time, reported:** 449 s and 379 s against 312 s and 282 s, at the flags' threads (the current method 8, each
  integrated step 2).
- **O2 is two-thirds of integrated CPU:** 533 s and 491 s.
- On round 16's cells the integrated stages are cheaper and faster: 119-168 CPU seconds against 211-287, and 68-100 s wall
  against 118-138.


---

## PREREG_o3_graph_round17_2026-10-17.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round17_2026-10-17.md`)*

## Round 17: integrated v13 on name-free families (linked by sequence identity and gene structure), three seeds per profile, and the cost of full families in CPU seconds (2026-10-17)

Pre-registration, written before any world of round 17's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 16 failed only on B1** (register 1318): in 3 of 36 worlds, 4-7 correct reads short of the current method, with no
  error where the current method has none.
- **v13 fixes the mechanisms on 52 spent runs** (register 1319), round 16's included:
  - thin and tied path calls are decided again on structure-aware distances;
  - reads with no path alignment are realigned.
- **The families are name-free** (user, 2026-10-17: "we should go by families linked by sequence identity, also if possible
  by intron/exon structure"). The advisor does not want the annotation leaned on, and the name-stem pools of rounds 11-16
  missed unnamed paralogs.

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v13** = v12 (register 1317) + two O2 genome-step rules (68722272, 67ba1cd2). On 52 spent runs: 44 -> 46 runs pass,
no world newly fails, correct +4,289, errors +95 (OUTSIDE 84, BRANCH 11).
- **`--struct-rescue`:** a thin or tied path call is decided again on gapstruct path distances (a long indel at its length).
- **`--short-path-align`:** reads the path step left without any alignment are aligned again with minimap2 -s 40.

**Integrated flags (`$V13` below):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

### Families (held out; name-free; drawn by a rule fixed before the build)

**The pool** (`seqfamilies.py`, 9b985ec5; structure test amended before the draw):
- **Genes:** every CAT protein-coding gene, whatever its name, with span <= 60 kb and overlaps resolved: 12,583 genes.
- **Sequence link:** gene A's canonical transcript, spliced-aligned to gene B's locus (+- 1 kb), covers >= 0.8 of A at identity
  >= 0.90, in either direction (families.py's locus rule, genome-wide). That gives 4,738 links among 1,119 genes.
- **Structure link:** one gene's exon chain is carried by the other's genome. In that spliced alignment every intron lies within
  10 bp of one of A's junctions, and none of A's junctions is left without one, in either direction. 4,384 of the 4,738 links
  pass; known retrocopies fail (SETSIP and SET, NUDT4B and NUDT4).
- **The amendment, before the draw:** the first version compared A's alignment with B's annotated introns. That failed identical
  genes whose canonical models differ (OR4F21 and OR4F16; ATAD3A's extra exon): a choice of the annotation, not of the genomes.
- **Families:** connected components of the links that pass both tests, 244 of them; 128 are unspent (no gene in any earlier
  simulated set), 282 copies.

**The draw** (fixed in 9b985ec5, before the build): the unspent families sorted by name, `random.Random(17).shuffle`, taken in
that order until the copy total reaches 90, a family skipped if it would take the total past 150.
- **40 families, 91 copies:** 33 pairs, 6 trios and one family of 7.
- **Four are X/Y pseudoautosomal pairs** (ASMTL, PPP2R3B, SHOX, SLC25A6), identical by construction.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families17.json` (set `r17`, sha1 b869562e5628), and `r17seq/families17seq_all.json`
  (+ `.links.tsv`).

| family | copies | genes | chromosomes |
|---|---|---|---|
| S17ADH1A | 2 | ADH1A, ADH1B | chr4 |
| S17AL6272301 | 7 | AL391987.2, AL627230.1 | chr9, chr13, chr21 |
| S17ASMTL | 2 | ASMTL | chrX, chrY |
| S17CBWD2 | 2 | CBWD1, CBWD2 | chr2, chr9 |
| S17CXorf49 | 2 | CXorf49, CXorf49B | chrX |
| S17DEFB104A | 3 | DEFB104A, DEFB104B | chr8 |
| S17DEFB107A | 3 | DEFB107A, DEFB107B | chr8 |
| S17DEFB4B | 3 | DEFB4A, DEFB4B | chr8 |
| S17DHRS4 | 2 | DHRS4, DHRS4L2 | chr14 |
| S17FCGR3A | 3 | FCGR3A, FCGR3B | chr1 |
| S17H2BC21 | 2 | H2AC20, H2BC21 | chr1 |
| S17H4C14 | 2 | H4C14, H4C15 | chr1 |
| S17KRT6B | 3 | KRT6A, KRT6B, KRT6C | chr12 |
| S17LCE3E | 2 | LCE3D, LCE3E | chr1 |
| S17LGALS7 | 2 | LGALS7, LGALS7B | chr19 |
| S17LIMS3 | 2 | LIMS3, LIMS4 | chr2 |
| S17LYZL1 | 2 | LYZL1, LYZL2 | chr10 |
| S17METTL2A | 2 | METTL2A, METTL2B | chr7, chr17 |
| S17MT1E | 2 | MT1E, MT1M | chr16 |
| S17NBPF6 | 2 | NBPF4, NBPF6 | chr1 |
| S17OR1L4 | 2 | OR1L4, OR1L6 | chr9 |
| S17OR1S2 | 2 | OR1S1, OR1S2 | chr11 |
| S17OR2T3 | 2 | OR2T3, OR2T34 | chr1 |
| S17OR52E6 | 2 | OR52E6, OR52E8 | chr11 |
| S17OR52K2 | 2 | OR52K1, OR52K2 | chr11 |
| S17OR5H1 | 2 | OR5H1, OR5H15 | chr3 |
| S17OR6B2 | 2 | OR6B2, OR6B3 | chr2 |
| S17OR9G1 | 2 | OR9G1 | chr11 |
| S17PPP2R3B | 2 | PPP2R3B | chrX, chrY |
| S17PRAMEF7 | 2 | PRAMEF7, PRAMEF8 | chr1 |
| S17RFPL4A | 2 | RFPL4A, RFPL4AL1 | chr19 |
| S17SFTPA2 | 2 | SFTPA1, SFTPA2 | chr10 |
| S17SHOX | 2 | SHOX | chrX, chrY |
| S17SLC25A6 | 2 | SLC25A6 | chrX, chrY |
| S17SMIM10L2B | 2 | SMIM10L2A, SMIM10L2B | chrX |
| S17SPACA5 | 2 | SPACA5, SPACA5B | chrX |
| S17SPATA31D4 | 3 | SPATA31D1, SPATA31D3, SPATA31D4 | chr9 |
| S17SPRR1A | 2 | SPRR1A, SPRR1B | chr1 |
| S17ULBP2 | 2 | RAET1L, ULBP2 | chr6 |
| S17XCL2 | 2 | XCL1, XCL2 | chr1 |

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`.
- **Seeds:** 701, 702 and 703, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families17.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout17
for S in 701 702 703; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r17 $S --local $O/local $V13 --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

### Bars

**Correctness (as rounds 15-16).** Read outcomes are per-read detectability less each world's error floor, scored on the veto
table (`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 17 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as round 16.** Measured on the mainstay set (NPIP as one family, TBC1D3, GOLGA8, GOLGA6L, CT47A,
XAGE, SPATA31; seed 306; both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r17cost/v13/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V13 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done
```
- **Bar:** v13's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, and the CPU seconds of round 17's cells.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that v13 does at least as well as the current method in every world type, on 40 held-out
families chosen by sequence and structure rather than by name, under both profiles and three seeds.

**A cost PASS** would claim that, on NPIP-sized families, the integrated pipeline uses no more CPU than the current method.

**Neither would claim:**
- anything on real reads;
- that near-identical uncatalogued twins are handled (register 1314-1315);
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs.

The pipeline under test still starts from the annotation's copies and isoforms (guided O1); round 17 removes the gene names from
the choice of test families, not the annotation from the method.

### Results (2026-10-17 23:34 - 2026-10-18 00:14; nothing changed after the first run, no crash)

**Verdict: FAIL, 0 of 6 cells. C1 PASSES** in both profiles.

| cell | failing worlds (integrated vs current: correct / errors) |
|---|---|
| A119b 701 | OUTSIDE 3,727 / 668 vs 893 / 191 |
| public 701 | OUTSIDE 3,671 / 781 vs 931 / 173 |
| A119b 702 | OUTSIDE 3,690 / 367 vs 725 / 163 |
| public 702 | OUTSIDE 3,487 / 376 vs 769 / 120; REMOVED 2,988 / 2 vs 803 / 1 |
| A119b 703 | OUTSIDE 3,792 / 540 vs 531 / 109 |
| public 703 | OUTSIDE 3,641 / 549 vs 530 / 112 |

- **B1 holds in all 36 worlds.** Round 16's shortfall is gone: v13 makes 124,264 correct calls against the current method's
  82,008.
- **Every failing bar but one is OUTSIDE's B2.** Over the 36 worlds v13 makes fewer errors, 6,986 against 8,489, but in OUTSIDE
  it makes 3-5 times the current method's.
- **The mechanism, post hoc (A119b 701 OUTSIDE):** 639 of the 668 errors are plain path calls on a family whose uncatalogued
  twin O1 did not seed as a new copy, so the catalogue keeps only the other copy:
  - LGALS7: 223 errors; CXorf49: 197; SPACA5: 27; SLC25A6: 18. The current method makes none in these families.
  - The twin's reads go to the only catalogued copy as a sole candidate. The current method's genomic candidates tie the read
    with the twin's locus, and it abstains.
- **This is the O2 proof's open failure (register 1314-1315),** made dominant by a set of pairs: take one copy of a pair out of
  the catalogue, and its near-identical twin has no rival in O2.
- **The guard measured for it in register 1315** (abstain when a locus off every seeded copy explains the read as well) cut these
  errors to 0-2. It was rejected because it also abstained the reads of an expressed copy tied with a silent catalogued copy
  that O1 had dropped.

**C1, the cost on the mainstay set:** v13 used 830 CPU seconds (A119b) and 675 (public) against the current method's 1,042 and
952 (-20% and -29%). Wall time, reported: 469 s and 378 s against 313 s and 281 s, at the flags' threads. Round 17's cells
used 591 CPU seconds in all against 1,033.


---

## PREREG_o3_graph_round18_2026-10-18.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round18_2026-10-18.md`)*

## Round 18: integrated v14 (copies the catalogue lacks, found in the genome) on fresh name-free families, three seeds per profile, and the cost in CPU seconds (2026-10-18)

Pre-registration, written before any world of round 18's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 17 failed every cell on OUTSIDE's B2** (register 1320). When a pair's twin is uncatalogued, O1 seeds nothing for it
  (no read is uniquely placed), and O2 gave the twin's reads to the one remaining copy.
- **The user chose to make it O3's job** (2026-10-18). Before O2, the catalogue is completed from the genome: a locus that holds a
  catalogued copy by the family rule becomes a new copy.
- **On round 17's six cells, now spent:** 0 -> 5 cells pass, and OUTSIDE errors fall from 367-781 to 1-159 (the current method:
  109-191).
  - B1 holds in all 36 worlds, at 1-8 correct reads a world.
  - The sixth cell fails REMOVED by one read (2 errors against 1), as v13 did.

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v14** = v13 (register 1319) + O1 `--ref-copies` (24bd5cde, 581570cf).
- Each kept curated copy's representative transcript is spliced-aligned to the world's homology-closed local reference (the
  masked one in REMOVED).
- An off-catalogue locus that holds it by the family rule becomes a new copy: identity >= 0.90 over >= 0.8 of it, and its exon
  chain carried by the locus, a junction at an unaligned end needing no intron.
- These copies get no `--twins` annotation prior. A read-seeded new copy that the search confirms loses that prior too.

**Integrated flags (`$V14`):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

### Families (held out; name-free; drawn by a rule fixed before the pool was rebuilt)

**The pool:** `seqfamilies.py build` (51f9b3a6), round 17's rule with the amended structure test (581570cf: a junction at an
unaligned end needs no intron).
- 4,508 of 4,738 sequence links pass structure, forming 256 families.
- **Spent** are all families touching a gene of any earlier simulated set, round 17's included. That leaves 98 unspent families,
  210 copies.

**The draw** (fixed in 51f9b3a6, before the rebuild): the unspent families sorted by name, `random.Random(18).shuffle`, taken in
that order until the copy total reaches 90, a family skipped if it would take the total past 150.
- **43 families, 91 copies;** copies per family: {2: 39, 3: 3, 4: 1}.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families18.json` (set `r18`, sha1 e689764bcf06), and `r18seq/families18seq_all.json`.

| family | copies | genes | chromosomes |
|---|---|---|---|
| S17ACSM2A | 2 | ACSM2A, ACSM2B | chr16 |
| S17AKAP17A | 2 | AKAP17A | chrX, chrY |
| S17AKR1B10 | 2 | AKR1B10, AKR1B15 | chr7 |
| S17AKR1C1 | 2 | AKR1C1, AKR1C2 | chr10 |
| S17AKR7A3 | 2 | AKR7A2, AKR7A3 | chr1 |
| S17ALG10 | 2 | ALG10, ALG10B | chr12 |
| S17ASMT | 2 | ASMT | chrX, chrY |
| S17ATAD3B | 2 | ATAD3A, ATAD3B | chr1 |
| S17CCDC74B | 2 | CCDC74A, CCDC74B | chr2 |
| S17CCL4L2 | 2 | CCL4L2 | chr17 |
| S17CD99 | 2 | CD99 | chrX, chrY |
| S17CELA3B | 2 | CELA3A, CELA3B | chr1 |
| S17CRLF2 | 2 | CRLF2 | chrX, chrY |
| S17CYP4A11 | 2 | CYP4A11, CYP4A22 | chr1 |
| S17EIF3C | 2 | EIF3C, EIF3CL | chr16 |
| S17ERVV1 | 2 | ERVV-1, ERVV-2 | chr19 |
| S17FAM86B2 | 2 | FAM86B1, FAM86B2 | chr8 |
| S17FFAR3 | 2 | FFAR3, GPR42 | chr19 |
| S17GTF2IRD2 | 2 | GTF2IRD2, GTF2IRD2B | chr7 |
| S17HTN3 | 2 | HTN1, HTN3 | chr4 |
| S17IFITM2 | 2 | IFITM2, IFITM3 | chr11 |
| S17IL9R | 2 | IL9R | chrX, chrY |
| S17MAGEB6B | 2 | MAGEB6, MAGEB6B | chrX |
| S17MAGED4B | 2 | MAGED4, MAGED4B | chrX |
| S17MTRNR2L6 | 2 | MTRNR2L6 | chr7, chr9 |
| S17OCM | 2 | OCM, OCM2 | chr7 |
| S17OR10H2 | 3 | OR10H1, OR10H2, OR10H5 | chr19 |
| S17OR4F21 | 2 | OR4F16, OR4F21 | chr1, chr5 |
| S17OR52I2 | 2 | OR52I1, OR52I2 | chr11 |
| S17OR8B2 | 2 | OR8B2, OR8B3 | chr11 |
| S17PLGLB1 | 3 | PLG, PLGLB1, PLGLB2 | chr2, chr6 |
| S17PRAMEF1 | 4 | PRAMEF1, PRAMEF13, PRAMEF14, PRAMEF2 | chr1 |
| S17PRAMEF18 | 2 | PRAMEF18, PRAMEF19 | chr1 |
| S17RAB6C | 2 | RAB6C, RAB6D | chr2 |
| S17RFPL1 | 3 | RFPL1, RFPL2, RFPL3 | chr22 |
| S17ROPN1 | 2 | ROPN1, ROPN1B | chr3 |
| S17SYT15 | 2 | FO681492.1, SYT15 | chr10 |
| S17TCP11X1 | 2 | TCP11X1, TCP11X2 | chrX |
| S17TRIM16 | 2 | TRIM16, TRIM16L | chr17 |
| S17TRIM43B | 2 | TRIM43, TRIM43B | chr2 |
| S17TUBB2A | 2 | TUBB2A, TUBB2B | chr6 |
| S17USP18 | 2 | USP18, USP41 | chr22 |
| S17ZBED1 | 2 | ZBED1 | chrX, chrY |

### Worlds and runs

- **Worlds:** six world types under `profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json`.
- **Seeds:** 801, 802 and 803, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families18.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout18
for S in 801 802 803; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r18 $S --local $O/local $V14 --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

### Bars

**Correctness (as rounds 15-17).** Per-read detectability less each world's error floor, on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 18 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as rounds 16-17,** on the mainstay set (seed 306, both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r18cost/v14/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V14 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done
```
- **Bar:** v14's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, the cells' CPU seconds, and the genome copies O1 found per world.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that v14 does at least as well as the current method in every world type, on 43
held-out name-free families, under both profiles and three seeds.

**A cost PASS** would claim no more CPU than the current method on NPIP-sized families.

**Neither would claim:**
- anything on real reads;
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs.

The catalogue still starts from the annotation's copies; the genome search adds the copies it lacks.

### Results (2026-10-18 10:57-11:43; nothing changed after the first run, no crash)

**Verdict: FAIL -- 4 of 6 cells pass; A119b 803 and public 803 fail. C1 PASSES** in both profiles.

| cell | failing worlds (integrated vs current: correct / errors) |
|---|---|
| A119b 803 | CROSS 4,052 / 0 vs 4,093 / 33; INTACT 3,726 / 0 vs 3,771 / 34 |
| public 803 | CROSS 3,703 / 0 vs 3,799 / 42; INTACT 3,391 / 0 vs 3,498 / 29; SILENT 2,140 / 0 vs 2,226 / 21 |

- **B2 holds in all 36 worlds.** Over the 36 worlds v14 makes 130,722 correct calls against the current method's 86,384, and
  3,782 errors against 10,431.
- **Round 17's failure is gone on fresh families:** OUTSIDE errors are 0-11 against the current method's 69-198, and REMOVED
  errors 1-17 against 26-85.
- **Every failing bar is B1, all in seed 803.** In those worlds v14 makes no errors and the current method makes 21-42.

**The mechanism, post hoc: one family, S17OR4F21 (OR4F21 on chr1, OR4F16 on chr5).**
- **O1's genome search adds five copies to it:** OR4F2P (chr11), OR4F7P (chr6), OR4F8P (chr19), OR4F28P (chr15) and
  FO082796.2 (chr3), 98.8-99.4% identical to OR4F16's transcript. All are unprocessed pseudogenes in CAT, and they are 5 of the 13
  genome copies O1 finds in each INTACT world.
- **The simulator expresses only the catalogue's protein-coding copies,** so these five never have a read.
- **Most reads the current method gets right and v14 does not are this family's:** 102 of 111 in A119b 803 INTACT, 159 of 172 in
  public 803 INTACT.
  - Tied with a pseudogene copy: 64 and 101.
  - Inside the fixed 2-edit margin of one (`thin`, mostly one substitution): 32 and 57.
- **The current method runs with `--no-as-tied-only`** (O3 needs it). That turns off the shipped tie-outside-catalogue rule
  (2026-09-09), so copy_assign compares only the catalogued copies.
  - Read s0000757's primary alignment is on OR4F2P, at the same AS as on OR4F21 (424).
  - copy_assign assigns it to OR4F21 at p 1e-18, against OR4F16 alone.
- **Seed 803 drew OR4F21 at 253-286 reads a world,** against 28-29 in seeds 801-802, where the same mechanism costs a few reads
  (B1 margins of +4 to +20 in INTACT).
- **In the real A119b testis library** (CHM13 v2.0 primary alignments), OR4F21 has 33 reads and OR4F16 26.
  - The five pseudogenes have 41: OR4F7P 18 (17 at MAPQ >= 1), OR4F8P 7, FO082796.2 7, OR4F28P 6, OR4F2P 3.
  - The simulation's silent pseudogenes are not what real testis shows.

**Two post-hoc recounts on these runs** (spent; not a verdict):
- **B1 on decidable reads only** (a catalogued read is undecidable when one of its alignments elsewhere in E0.bam scores at least
  its best at its true copy): 4 of the 5 failing bars pass. Public 803 INTACT stays short, 3,383 against 3,392, by the
  one-substitution reads.
- **Thin reads decided by the calibrated significance gate alone,** in place of the fixed 2-edit margin (offline replay, veto and
  structure check not replayed):
  - +100 to +566 correct a world, with no added error in INTACT or SILENT and 2-7 in CROSS.
  - OUTSIDE is not scored by the replay, which does not map new copies to their truth copies.
  - Its genome-step form, `--significance-margin`, was rejected in register 1319: one error per world where the current method
    makes none.

**C1, the cost on the mainstay set:** v14 used 1,014 CPU seconds (A119b) and 820 (public), against the current method's 1,044 and
945 (-3% and -13%).
- O1's genome search costs 141-147 CPU seconds of these.
- Reported: wall 565 s and 473 s against 310 s and 283 s, at the flags' threads; peak RSS 2.42 and 2.35 GB.
- Round 18's cells used 1,336 CPU seconds in all against 1,254. The set's one-time local-reference build peaks at 14.4 GB.


---

## PREREG_o3_graph_round19_2026-10-18.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round19_2026-10-18.md`)*

## Round 19: integrated v14 (unchanged) on fresh families whose non-coding members are copies like any other, expressed as real testis expresses them (2026-10-18)

Pre-registration, written before any world of round 19's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 18 failed the two seed-803 cells on B1** (register 1322; 4 of 6 cells passed). Nearly every lost read was OR4F21's: O1's
  genome search found five OR4F unprocessed pseudogenes, the simulator never expressed them, and v14 abstained on the reads
  tied with them. The current method compares catalogued copies only, so in that simulation it was always right.
- **Real A119b testis expresses those pseudogenes** (41 primary reads against 59 on the two coding copies).
- **The user's decision (2026-10-18):** non-coding paralogs are family members like any other. They go into the catalogue for
  both methods, and the simulator expresses them at the level the real libraries show. Retrocopies the family rule admits (those
  of single-exon genes, where structure cannot decide) are kept.
- **Development on round 18's families with their non-coding members** (register 1323, spent): 6 of 6 cells pass. The thinnest
  B1 margin is +93 reads (round 18: -107) and the thinnest B2 margin +7.

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v14:** unchanged since round 18 (register 1321), with the same flags:
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

**Only the simulation changed:**
- **Families:** `seqfamilies.py build --paralogs` (d8ea648d). A gene of any biotype joins a family when its locus carries a coding
  gene's canonical transcript by the family rule (identity >= 0.90 over >= 0.8, exon chain carried; that direction only).
  - A retrocopy of a multi-exon gene fails the structure test.
  - None may overlap a coding gene, and a pseudogene wins over a lncRNA that overlaps it.
- **Expression:** profiles v4 (19a1b482) = v2 plus the measured expression of non-coding members (`calibrate.py paralogs`,
  uniquely placed reads, families whose coding copies average >= 10 reads):
  - A119b: 3.7% silent, the rest expressed at a median 0.34 of the family's coding copies;
  - public testis: 37.5% silent, the rest at 0.14.
  - A silent member stays catalogued (kind `silent`).
- **Scenario copies** (REMOVED, OUTSIDE, BRANCH) and SILENT's silenced third are drawn among the coding copies, as before
  (d8ea648d; BRANCH's extra copy fixed in b2502ca4).

### Families (held out; drawn by the fixed rule)

**The pool:** `seqfamilies.py build --paralogs` with every earlier simulated set spent, round 18's included (23 files,
`ls families*.json | grep -v _all.json`).
- 341 families by sequence and structure; 130 are unspent (301 copies, 106 of them non-coding).
- **File:** `r19seq/families19seq_all.json`.

**The draw:** the rule of rounds 17-18 (51f9b3a6) with the round's number as seed: the unspent families sorted by name,
`random.Random(19).shuffle`, taken in that order until the copy total reaches 90, a family skipped if it would pass 150.
- **35 families, 90 copies** (54 coding, 36 non-coding); copies per family {2: 23, 3: 7, 4: 3, 5: 1, 6: 1}.
- **Non-coding members:** 19 unprocessed pseudogenes, 4 transcribed unprocessed, 8 processed, 3 transcribed processed, 1
  polymorphic, 1 unitary.
  - 9 of the 11 retrocopies are in families of single-exon genes.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families19.json` (set `r19`, sha1 73ee6bdb294c).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19ADAM21 | 2 | ADAM21 | ADAM21P1 (tr. proc.) | chr14 |
| S19ANXA8L1 | 2 | ANXA8L1, ANXA8 | - | chr10 |
| S19AZGP1 | 2 | AZGP1 | AZGP1P1 (tr. unproc.) | chr7 |
| S19DDX11 | 3 | DDX11 | AC009533.1 (unproc.), DDX12P (unproc.) | chr12 |
| S19DEFB106A | 3 | DEFB106A, DEFB106A, DEFB106B | - | chr8 |
| S19DEFB130A | 3 | DEFB130A | DEFB130C (tr. unproc.), DEFB130D (unproc.) | chr4, chr8, chr11 |
| S19GPAT2 | 2 | GPAT2 | GPAT2P1 (tr. unproc.) | chr2 |
| S19GPR32 | 2 | GPR32 | GPR32P1 (unproc.) | chr19 |
| S19H2BC9 | 4 | H2BC9, H2BC13, H2AC13, H2BC17 | - | chr6 |
| S19H4C11 | 2 | H4C11, H4C12 | - | chr6 |
| S19HSD3B2 | 4 | HSD3B2, HSD3B1 | HSD3BP2 (unproc.), HSD3BP4 (unproc.) | chr1 |
| S19IKBKG | 2 | IKBKG | IKBKGP1 (unproc.) | chrX |
| S19IL3RA | 2 | IL3RA, IL3RA | - | chrX, chrY |
| S19NACA2 | 3 | NACA2 | NACA3P (proc.), NACA4P (tr. proc.) | chr4, chr8, chr17 |
| S19NAT8 | 2 | NAT8 | NAT8B (polymorphic) | chr2 |
| S19NUDT4B | 2 | NUDT4B | NUDT4P2 (proc.) | chr1 |
| S19OR4A5 | 2 | OR4A5 | OR4A42P (unproc.) | chr11 |
| S19OR4K15 | 3 | OR4K15 | OR4K8P (unproc.), OR4K11P (unproc.) | chr14, chr18, chr21 |
| S19OR5H14 | 3 | OR5H14, OR5H6 | OR5H7P (unproc.) | chr3 |
| S19OR5M8 | 2 | OR5M8 | OR5M6P (unproc.) | chr11 |
| S19OR8I2 | 2 | OR8I2 | OR8I1P (unproc.) | chr11 |
| S19PABPC1L2B | 2 | PABPC1L2B, PABPC1L2A | - | chrX |
| S19PDZK1 | 2 | PDZK1 | PDZK1P1 (tr. unproc.) | chr1 |
| S19PGAM4 | 4 | PGAM4 | PGAM1P8 (tr. proc.), AC012085.1 (proc.), AC078819.1 (proc.) | chrX, chr11, chr12 |
| S19PRAMEF10 | 5 | PRAMEF10, PRAMEF33, PRAMEF17 | PRAMEF31P (unproc.), PRAMEF31P (unproc.) | chr1 |
| S19PRR23D1 | 6 | PRR23D1, PRR23D2, PRR23D1 | AC134684.4 (unproc.), PRR23D3P (unproc.), AC134684.4 (unproc.) | chr8 |
| S19PTTG2 | 2 | PTTG2 | PTTG3P (proc.) | chr4, chr8 |
| S19RHOXF2B | 2 | RHOXF2B, RHOXF2 | - | chrX |
| S19SAA2 | 2 | SAA2, SAA1 | - | chr11 |
| S19SERF1B | 2 | SERF1B, SERF1A | - | chr5 |
| S19SLC25A52 | 3 | SLC25A52 | SLC25A51P1 (proc.), SLC25A51P3 (proc.) | chr6, chr8, chr18 |
| S19TGIF2LX | 2 | TGIF2LX, TGIF2LY | - | chrX, chrY |
| S19VN1R4 | 2 | VN1R4 | VN1R6P (unitary) | chr19 |
| S19ZC3H11B | 2 | ZC3H11B | AL365217.1 (proc.) | chr1, chr6 |
| S19ZNF181 | 2 | ZNF181 | ZNF807P (unproc.) | chr19 |

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 901, 902 and 903, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families19.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout19
for S in 901 902 903; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r19 $S --local $O/local $V14 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

### Bars

**Correctness (as rounds 15-18).** Per-read detectability less each world's error floor, on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 19 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as rounds 16-18,** on the mainstay set (seed 306, both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r19cost/v14/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V14 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done
```
- The mainstay families carry no biotypes, so v4 simulates them as v2 did.
- **Bar:** v14's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, the cells' CPU seconds, and the genome copies O1 found per world.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that v14 does at least as well as the current method in every world type, on 35 held-out
name-free families with their non-coding members expressed at the measured level, under both profiles and three seeds.

**A cost PASS** would claim no more CPU than the current method on NPIP-sized families.

**Neither would claim:**
- anything on real reads;
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs;
- that the family definition is settled: the Compara study (in design, 2026-10-18) may change it;
- anything about paralogs no annotation holds: loci found only by O1's genome search still express nothing.

### Results (2026-10-18 16:13-16:52; nothing changed after the first run, no crash)

**Launch note:** the first launch (16:12) failed before any run started. The batch's log directories (`fast/r19cost/`,
`fast/heldout19/`) did not exist, so every run's output redirect failed and nothing was simulated. The directories were
created and the same batch relaunched unchanged.

**Verdict: PASS -- 6 of 6 cells, all 72 bar checks. C1 PASSES** in both profiles.

| cell | INTACT | CROSS | OUTSIDE | REMOVED | BRANCH | SILENT |
|---|---|---|---|---|---|---|
| A119b 901 | 2,785 / 0 vs 2,682 / 0 | 2,828 / 0 vs 2,725 / 0 | 2,340 / 1 vs 722 / 780 | 2,375 / 0 vs 1,355 / 24 | 3,664 / 516 vs 2,743 / 1,035 | 1,461 / 0 vs 1,402 / 0 |
| public 901 | 2,359 / 0 vs 2,213 / 21 | 2,382 / 0 vs 2,237 / 21 | 1,856 / 0 vs 433 / 853 | 1,827 / 0 vs 1,116 / 7 | 3,145 / 631 vs 2,149 / 1,212 | 1,063 / 0 vs 974 / 21 |
| A119b 902 | 2,725 / 0 vs 2,652 / 0 | 2,778 / 0 vs 2,705 / 0 | 2,342 / 0 vs 726 / 534 | 2,391 / 0 vs 1,144 / 25 | 3,295 / 256 vs 2,742 / 491 | 1,336 / 0 vs 1,285 / 0 |
| public 902 | 2,275 / 0 vs 2,086 / 0 | 2,294 / 0 vs 2,105 / 0 | 1,636 / 1 vs 461 / 587 | 1,856 / 1 vs 890 / 11 | 2,746 / 267 vs 2,174 / 424 | 874 / 0 vs 836 / 0 |
| A119b 903 | 3,149 / 0 vs 3,025 / 0 | 3,213 / 0 vs 3,089 / 0 | 2,846 / 0 vs 724 / 336 | 2,758 / 1 vs 1,043 / 8 | 3,918 / 324 vs 3,025 / 694 | 1,159 / 0 vs 1,121 / 0 |
| public 903 | 2,703 / 0 vs 2,530 / 0 | 2,739 / 0 vs 2,566 / 0 | 1,927 / 6 vs 434 / 295 | 2,261 / 1 vs 715 / 11 | 3,417 / 350 vs 2,553 / 638 | 717 / 0 vs 688 / 0 |

(integrated correct / errors vs current correct / errors)

- **Over the 36 worlds:** v14 makes 85,440 correct calls against the current method's 62,070, and 2,355 errors against 8,028.
- **B1's thinnest margin is +29 reads** (public 903 SILENT), then +38 (A119b 903 SILENT, public 902 SILENT).
- **B2's thinnest margins are ties at zero:** in 15 complete-catalogue worlds (INTACT, CROSS, SILENT) both methods make no error.
- **Genome copies O1 found per world:** 4-5 in INTACT, CROSS, REMOVED, BRANCH and SILENT, and 31-34 in OUTSIDE (seed 901).

**C1, the cost on the mainstay set:** v14 used 995 CPU seconds (A119b) and 825 (public), against the current method's 1,032 and
943 (-4% and -13%).
- Reported: round 19's cells used 849 CPU seconds in all against 750. On these small families O1's genome search (about 20
  CPU seconds a world) costs more than it saves.
- The set's one-time local-reference build peaks at 14.4 GB.

**What it claims, as registered:** v14 does at least as well as the current method in every world type, on 35 held-out name-free
families with their non-coding members expressed at the measured level, under both profiles and three seeds, and uses no more
CPU on NPIP-sized families. It claims nothing on real reads, on arrays of FAM90A's or DUX4's size, or on paralogs no annotation
holds.


---

## PREREG_o3_graph_round20_2026-10-18.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round20_2026-10-18.md`)*

## Round 20: integrated v15 (O1's genome search under the carried intron chain, every chain kept) on fresh families built by the same rule (2026-10-18)

Pre-registration, written before any world of round 20's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 19 passed** (register 1324): v14 held every bar on fresh families whose non-coding members are copies.
- **The family rule missed known families** (register 1325):
  - GOLGA6L failed the exact structure test, because copies split a long exon with short introns of their own and indels come
    out of the spliced alignment as introns.
  - TSPY2 was lost to minimap2's default chain selection among about 45 near-identical targets.
- **The user's decisions (2026-10-18):**
  - build the next family pool with the carried intron chain (`--intron-chain carried`) and every chain kept
    (`--retain-chains`);
  - apply both to O1's genome search.
- **Development of v15 on spent runs** (register 1326):
  - Round 19's six cells and the mainstay set, both profiles: 8 of 8 pass.
  - Cold-cache cost with the search cached per set (O1_ALN_CACHE): 1,019 / 848 CPU s against the current method's 1,031 / 955.
    A119b's margin is only 12 s.

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v15** = v14 + O1 `--ref-carried` (05cb9c6d, cached through O1_ALN_CACHE under `--aln-cache`, 5b3d50ac). O1's genome
search takes a locus as a copy when the catalogued transcript's intron chain is carried there (each junction has its intron; the
locus may hold more), with every minimap2 chain kept (-P).

**Integrated flags (`$V15`):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --ref-carried --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

**The simulation, as round 19:**
- profiles v4: non-coding members silent at the measured fraction, else expressed at the measured ratio;
- scenario copies drawn among the coding copies.
- **Fixed since round 19:** a world's regions end at their contig's end (49323d8b). Round 20's draw holds two mitochondrial
  families (chrM is 16,569 bp). No region of rounds 15-19 ran past an end.

### Families (held out; drawn by the fixed rule)

**The pool:** `seqfamilies.py build --paralogs --intron-chain carried --retain-chains`.
- Every earlier simulated set is spent, round 19's included (24 files, `ls families*.json | grep -v _all.json`).
- 341 families by sequence and structure; 96 are unspent (218 copies, 73 of them non-coding).
- **File:** `r20seq/families20seq_all.json`.

**The draw:** the rule of rounds 17-19 (51f9b3a6) with the round's number as seed: the unspent families sorted by name,
`random.Random(20).shuffle`, taken in that order until the copy total reaches 90, a family skipped if it would pass 150.
- **38 families, 90 copies** (60 coding, 30 non-coding); copies per family {2: 30, 3: 5, 4: 1, 5: 1, 6: 1}.
- **Non-coding members:** 11 unprocessed pseudogenes, 6 transcribed unprocessed, 10 processed, 1 polymorphic, 2 StringTie
  models.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families20.json` (set `r20`, sha1 fcf105530369).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19BX6646152 | 2 | BX664615.2 | FGF7P6 (tr. unproc.) | chr9 |
| S19CFC1B | 2 | CFC1B, CFC1 | - | chr2 |
| S19FAHD2A | 3 | FAHD2A, FAHD2B | FAHD2CP (tr. unproc.) | chr2 |
| S19FKBP9 | 2 | FKBP9 | FKBP9P1 (tr. unproc.) | chr7 |
| S19GK3P | 2 | GK3P | AC110373.1 (unproc.) | chr4 |
| S19GLYATL1B | 2 | GLYATL1B | GLYATL1P4 (unproc.) | chr11 |
| S19HCAR2 | 2 | HCAR2, HCAR3 | - | chr12 |
| S19HNRNPA1L2 | 4 | HNRNPA1L2 | HNRNPA1P7 (proc.), HNRNPA1P10 (proc.), HNRNPA1P35 (proc.) | chr2, chr13, chr18, chr19 |
| S19HP | 2 | HP, HPR | - | chr16 |
| S19HPF1 | 2 | HPF1 | AL592293.1 (unproc.) | chr4, chr9 |
| S19HSD17B7 | 2 | HSD17B7 | HSD17B7P2 (tr. unproc.) | chr1, chr10 |
| S19HSFX2 | 2 | HSFX2, HSFX1 | - | chrX |
| S19IFNA13 | 2 | IFNA13, IFNA1 | - | chr9 |
| S19IFNL3 | 2 | IFNL3, IFNL2 | - | chr19 |
| S19KRT81 | 2 | KRT81, KRT83 | - | chr12 |
| S19LAIR1 | 2 | LAIR1, LAIR2 | - | chr19 |
| S19MTCO2 | 2 | MT-CO2 | MTCO2P22 (proc.) | chr5, chrM |
| S19MTND5 | 5 | MT-ND5, MT-ND6 | MTCYBP18 (proc.), MTND6P4 (proc.), MTND5P11 (proc.) | chr5, chrM |
| S19MYH1 | 2 | MYH1, MYH2 | - | chr17 |
| S19OR10A5 | 2 | OR10A5, OR10A2 | - | chr11 |
| S19OR10H3 | 2 | OR10H3, OR10H4 | - | chr19 |
| S19OR2T27 | 2 | OR2T27 | OR2T7 (polymorphic) | chr1 |
| S19OR4C12 | 3 | OR4C12 | OR4C2P (unproc.), OR4C50P (unproc.) | chr11 |
| S19OR4C3 | 2 | OR4C3 | OR4C49P (unproc.) | chr11 |
| S19OR5AK2 | 2 | OR5AK2 | OR5AK1P (unproc.) | chr11 |
| S19OR5L1 | 2 | OR5L1, OR5L2 | - | chr11 |
| S19PDPR | 3 | PDPR | MSTRG.2242 (StringTie), MSTRG.2243 (StringTie) | chr16 |
| S19PNMA6A | 2 | PNMA6A | PNMA6B (unproc.) | chrX |
| S19RAB40AL | 2 | RAB40AL, RAB40A | - | chrX |
| S19REG3G | 2 | REG3G, REG3A | - | chr2 |
| S19RNFT1 | 2 | RNFT1 | RNFT1P3 (unproc.) | chr17 |
| S19RPL10L | 2 | RPL10L | RPL10P11 (proc.) | chr7, chr14 |
| S19SDHA | 2 | SDHA | SDHAP2 (tr. unproc.) | chr3, chr5 |
| S19TRIM74 | 3 | TRIM74, TRIM50, TRIM73 | - | chr7 |
| S19TUBA1B | 2 | TUBA1B, TUBA1A | - | chr12 |
| S19UBTFL1 | 3 | UBTFL1 | UBTFL7 (proc.), UBTFL2 (proc.) | chr11 |
| S19ZNF479 | 6 | ZNF479, ZNF722P, ZNF735, ZNF679 | ZNF733P (tr. unproc.), ZNF734P (unproc.) | chr7 |
| S19ZNF736 | 2 | ZNF736 | AC095038.1 (unproc.) | chr7 |

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 1001, 1002 and 1003, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families20.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout20
for S in 1001 1002 1003; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r20 $S --local $O/local $V15 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

### Bars

**Correctness (as rounds 15-19).** Per-read detectability less each world's error floor, on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 20 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as rounds 16-19,** on the mainstay set (seed 306, both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r20cost/v15/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V15 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done
```
- **Bar:** v15's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, the cells' CPU seconds, and the genome copies O1 found per world.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that v15 does at least as well as the current method in every world type, on 38 held-out
families built by the carried-intron-chain rule with their non-coding members expressed at the measured level, under both
profiles and three seeds.

**A cost PASS** would claim no more CPU than the current method on NPIP-sized families.

**Neither would claim:**
- anything on real reads;
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs;
- anything about paralogs no annotation holds.

### Results (2026-10-18 18:29-19:18; nothing changed after the first run, no crash)

**Verdict: FAIL -- 5 of 6 cells pass; A119b 1001 fails. C1 PASSES** in both profiles.

| cell | failing worlds (integrated vs current: correct / errors) |
|---|---|
| A119b 1001 | INTACT 3,287 / 3 vs 3,285 / 0; CROSS 3,407 / 3 vs 3,401 / 1 |

- **Over the 36 worlds:** v15 makes 111,360 correct calls against the current method's 78,757, and 2,857 errors against 8,644.
- **Every failing bar is B2, from one family and one mechanism:** in each world, three reads of ZNF479 (copy 0) are moved to
  ZNF734P (copy 2, an unprocessed pseudogene) by the genome veto (`veto_locus`).
  - Read s0003469's primary alignment lies on ZNF734P at MAPQ 36 and AS 700.
  - Its secondary on its own copy scores AS 724.
  - minimap2 chose the primary by its DP score (768 against 748). Its two introns, 84 and 252 bp, skip one and three
    zinc-finger repeat units of 84 bp. The alignment at the read's own locus has 9 clustered mismatches in the repeats.
  - The veto trusts a primary at MAPQ >= 30. It does not ask whether a secondary on another copy scores as well, which the
    user's rule of 2026-10-11 forbids: a primary with an AS-tied secondary is not trusted as the origin.
- **The weakness predates v15:** under round 19's exact rule, ZNF479's family is exactly {ZNF479, ZNF734P}. Round 20 is the
  first draw to hold it.

**C1, the cost on the mainstay set:** v15 used 1,028 CPU seconds (A119b) and 854 (public), against the current method's 1,046
and 945.
- O1's genome search, cached per set, costs 170 s of these.
- Reported: round 20's cells used 1,311 CPU seconds in all against 1,078.


---

## PREREG_o3_graph_round21_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round21_2026-10-19.md`)*

## Round 21: O2 only on ambiguous reads, by Eichler's rule (within 10 AS), on fresh multigene families (2026-10-19)

This pre-registration was written before any world of round 21's families was built; no read of these families has been
simulated. Code: branch `o3-fast`, frozen at the commit that adds this file.

### Why

**The user (2026-10-19):**
- "go ahead and use only for tied AS reads"
- "only AS tied reads are ambiguous and need provenance; O2 on mapq>0 reads is tautological since they are already placed
  at their right placement so only multimappers with tied AS where the distinction between primary and secondary is
  meaningless need O2"
- "according to the Eichler paper, their definition of a well placed multimapper is that where the AS is within 10"

**Eichler's rule, as implemented here before** (register 935): a read is placed on its best-AS alignment iff no other
alignment scores within 10 AS units of it.
- Under splice:hq a substitution moves AS by 5 (+1 match, -4 mismatch), so 10 is two sequence differences.
- On HG002's untied reads (development, register 1330 data), minimap2's copy at a margin of 5-9 was wrong 46.5% of the time,
  and at a margin of 10 or more 1.6% of the time.

**Development (all spent):**
- **HG002 (real reads, the donor's own copies as truth):** this system and v14 agree.
  - On 49,057 untied reads with a donor copy: 44,064 right and 597 wrong, against v14's 44,021 and 599.
  - On the 582 tied reads: 49 right and 7 wrong, against 50 and 7.
  - Equal-AS-only scoping instead was wrong on 3.2%.
- **Round 19's 36 worlds, rescored with the scorer below:** 4 of 6 cells pass; v14 passes all 6.
  - Both failures are B2: 1-2 PABPC1L2B reads in each of three worlds. Their own copy carries the paralog's alleles, so
    minimap2 places them at a margin of 12 on the other copy, where v14's paths abstained.
  - The user kept B2 strict and accepted these reads as the split's known cost.
- **Two fixes found in development, both in the code:**
  - A catalogued copy O1 dropped as silent ties like a copy, not like a locus off the catalogue (register 1315's lesson).
  - Reads tied within 10 AS with a locus on no copy of the world are scored apart (user's decision): the simulation never
    expresses that locus.

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v16** = v14 (round 19, register 1324) + `--o2-scope tied --as-slack 9` (`o2graph.py run W intega`). Per read and
family, the alignments in `E0.bam` within 9 AS of the read's best are its loci:
- **No alignment:** abstain.
- **A locus off the family's copies:** abstain.
- **One seeded copy:** that copy, minimap2's placement. The exception is when its alignment does not explain the read (a
  novel insertion or clip, or divergence beyond DELTA plus the world's error floor); then abstain, leaving it to O3.
- **The silent copy alone:** abstain.
- **Two or more copies (a tie):** v14's call when it is one of them; otherwise abstain.

O3 is unchanged: it reads `o2.integ`, and its calls take precedence.

**Integrated flags (`$V16`):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
--o2-scope tied --as-slack 9
```

**The simulation:** as rounds 19-20 (profiles v4).

### Families (held out; drawn by the fixed rule)

**The pool:** the adopted multigene rule (`seqfamilies.py build --rule multigene`, register 1328), with all components
(`r19seq/families19cal_full2_all.json`, 439 families).
- **Spent:** every family touching a gene of an earlier simulated set (26 files, `ls families*.json | grep -v _all.json`),
  or of the HG002 worlds (`hg002proof/world*/families.json`, seen in development).
- That leaves 135 unspent families with 393 copies. File: `r21seq/families21seq_all.json`.

**The draw:** the rule of rounds 17-20 with seed 21: names sorted, `random.Random(21).shuffle`, taken in order until the
copy total reaches 90, skipping a family that would pass 150.
- **37 families, 90 copies:** 45 coding and 45 non-coding (33 transcribed unprocessed, 8 unprocessed, 3 transcribed
  processed, 1 processed).
- **Copies per family:** {2: 28, 3: 4, 4: 4, 6: 1}.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families21.json` (set `r21`, sha1 0e002f85636e).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19AC0928653 | 6 | - | AC092865.3 (tr. unproc.), AC079781.4 (tr. unproc.), AC068587.8 (tr. unproc.), AC068587.8 (tr. unproc.), AC105233.3 (tr. unproc.), AC068587.6 (tr. unproc.) | chr7, chr8, chr12 |
| S19AC1335553 | 3 | - | AC133555.3 (tr. unproc.), AC106782.1 (tr. unproc.), AC106782.1 (tr. unproc.) | chr16 |
| S19AC2455951 | 4 | - | AC245595.1 (tr. unproc.), AC245595.1 (tr. unproc.), AC027612.2 (tr. unproc.), BAGE2 (tr. unproc.) | chr1, chr2, chr21 |
| S19ADAMTS7P1 | 2 | - | ADAMTS7P1 (tr. unproc.), ADAMTS7P4 (tr. unproc.) | chr15 |
| S19ANAPC1 | 2 | ANAPC1 | ANAPC1P4 (tr. unproc.) | chr2 |
| S19AQP12B | 2 | AQP12B, AQP12A | - | chr2 |
| S19BMS1P8 | 4 | - | BMS1P8 (tr. unproc.), BMS1P8 (tr. unproc.), AC026273.1 (tr. unproc.), AC026273.1 (tr. unproc.) | chr2, chr16 |
| S19CCT8L2 | 2 | CCT8L2 | CCT8L1P (proc.) | chr7, chr22 |
| S19CD177 | 2 | CD177 | CD177P1 (unproc.) | chr19 |
| S19CSF2RA | 2 | CSF2RA | - | chrX, chrY |
| S19CYP21A2 | 2 | CYP21A2 | CYP21A1P (tr. unproc.) | chr6 |
| S19CYP2C19 | 2 | CYP2C19, CYP2C9 | - | chr10 |
| S19CYP3A7 | 2 | CYP3A7, CYP3A4 | - | chr7 |
| S19DHX40 | 2 | DHX40 | DHX40P1 (tr. unproc.) | chr17 |
| S19EOLA1 | 2 | EOLA1, EOLA2 | - | chrX |
| S19MALL | 2 | MALL | AC140479.5 (tr. unproc.) | chr2 |
| S19MSTRG2518 | 4 | MSTRG.2518 | USP32P1 (tr. unproc.), USP32P3 (tr. unproc.), AC053481.2 (tr. unproc.) | chr17 |
| S19MSTRG676 | 2 | MSTRG.676 | PARGP1 (tr. unproc.) | chr10 |
| S19NCOR1P4 | 3 | - | NCOR1P4 (tr. proc.), NCOR1P4 (tr. proc.), NCOR1P4 (tr. proc.) | chr13, chr14, chr21 |
| S19NLGN4X | 2 | NLGN4X, NLGN4Y | - | chrX, chrY |
| S19NXF2 | 2 | NXF2, NXF2B | - | chrX |
| S19OR1E1 | 2 | OR1E1, OR1E2 | - | chr17 |
| S19OR4A15 | 2 | OR4A15 | OR4A21P (unproc.) | chr11 |
| S19OR4C46 | 2 | OR4C46 | OR4C10P (unproc.) | chr11 |
| S19OR52L1 | 2 | OR52L1 | OR52L2P (unproc.) | chr11 |
| S19OR5F1 | 2 | OR5F1 | OR5F2P (unproc.) | chr11 |
| S19PDE4DIP | 3 | PDE4DIP | AC241952.1 (tr. unproc.), AC245297.1 (tr. unproc.) | chr1 |
| S19PROS1 | 2 | PROS1 | PROS2P (unproc.) | chr3 |
| S19RLN2 | 2 | RLN2, RLN1 | - | chr9 |
| S19S100A7A | 2 | S100A7A, S100A7 | - | chr1 |
| S19SORD | 2 | SORD | SORD2P (tr. unproc.) | chr15 |
| S19TBC1D26 | 2 | TBC1D26, TBC1D28 | - | chr17 |
| S19TEKT4P2 | 2 | - | TEKT4P2 (tr. unproc.), TEKT4P2 (tr. unproc.) | chr13, chr15 |
| S19TMEM191C | 2 | TMEM191C | TMEM191A (tr. unproc.) | chr22 |
| S19TPTE2 | 3 | TPTE2, TPTE | TPTE2P4 (unproc.) | chrY, chr13, chr21 |
| S19ZNF257 | 4 | ZNF257, ZNF724, ZNF675, ZNF273 | - | chr7, chr19 |
| S19ZNF285 | 2 | ZNF285 | AC018688.1 (unproc.) | chr19, chr20 |

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 1101, 1102 and 1103, with no seed selection. That makes 6 cells and 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families21.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout21
for S in 1101 1102 1103; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r21 $S --local $O/local $V16 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

### Bars

**Correctness.** As rounds 15-20: per-read detectability less each world's error floor, on the scoped table
(`--o2 o2.intega.m2.tsv`). Reads tied within 9 AS with a locus on no copy of the world (`--off-family-ties 9`) are the
class `off_family_tie`, left out of both systems' sums. fastset's own e2e step does this.
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**Verdict.** A cell PASSES when all six of its worlds pass both bars, and round 21 PASSES only when all six cells pass (72 bar
checks).

**Cost (C1), in CPU seconds,** as rounds 16-20, on the mainstay set (seed 306, both profiles), with nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r21cost/v16/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V16 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done
```
- **Bar:** v16's integrated CPU seconds (O1 + O2 + O2_genome + O2_scope + O3, summed over the six worlds) at most the
  current method's, in each profile.

**Reported, no bar:**
- v14 on the same worlds and scoring (`o2.integv.m2.tsv`, from the same runs).
- The `off_family_tie` reads' fates.
- The errors of reads minimap2 resolved (`as_resolved`) apart from O2's.
- Wall times and peak RSS.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that O2 confined to ambiguous reads (within 10 AS, Eichler's rule), with minimap2's
placement for the rest, does at least as well as the current method in every world type. It would hold on 37 held-out
multigene families, half of them with non-coding members.

**It would not claim:**
- **Anything about reads tied with a locus off every copy.** These are scored apart.
- **Accuracy where a donor's copy carries its paralog's alleles at two or more sites.** minimap2 places those reads on the
  paralog. Development shows such reads exist, and real data does too (register 1330).
- **Real-read accuracy beyond the HG002 development numbers above.**

### Results (2026-10-19, run once as pre-registered; `/mnt/linuxdisk/tmp/o3_graph/fast/heldout21/`)

The six cells ran in 62 min (7-9 min each) and the two cost runs in 22 min, with no crash. The code was unchanged from
097390a1.

**Verdict: FAIL.** 5 of 6 cells pass, and 69 of the 72 bar checks hold.
- **public_s1103 fails B2 in CROSS, INTACT and SILENT,** by one error each against the current method's 0.
- **The three errors are one read sequence:**
  - an NCOR1P4 read of copy .1 (chr14);
  - minimap2 aligns it to copy .2 (chr21) with AS 278 at NM 0, and to its own copy with AS 262, also at NM 0 (the primary,
    MAPQ 34);
  - correction (same day): the chr21 alignment covers LESS of the read. It leaves 38 bases soft-clipped and has 1 intron,
    while the own copy aligns all 340 bases across 3 introns. minimap2's AS charges each intron and nothing for a clip, so
    it favours the copy that explains less of the read. At a margin of 16 the read is minimap2's, and it goes to the wrong
    copy;
  - v14's paths placed it right.
- **This is the known weakness of AS margins** (register 1270): an alignment with fewer introns, its unexplained ends
  clipped for free, outscores the true copy's full spliced alignment. NCOR1P4's three copies are transcribed processed
  pseudogenes.

**C1 PASS** in both profiles. Integrated CPU, against the current method's:
- A119b: 938 s against 1,044. The O2 scope step took 5 s.
- public: 746 s against 952. The O2 scope step took 4 s.

**Reported, no bar:**
- **All 36 worlds, judged reads** (off-family ties left out):

  | system | correct | errors |
  |---|---|---|
  | v16 | 97,036 | 4,137 |
  | v14, same worlds and scoring | 97,132 | 4,162 |
  | current | 64,640 | 8,362 |

  v14 passes all six cells under this scoring.
- **Off-family ties** (scored apart):

  | system | correct | errors |
  |---|---|---|
  | v16 | 1,398 | 2 |
  | v14 | 1,746 | 37 |
  | current | 4,084 | 86 |

  These are the ties v14 and the current method guess on.
- **Who placed what** (INTACT, CROSS, OUTSIDE and SILENT; catalogued copies' reads; off-family ties left out):
  - minimap2 (`as_resolved`): 49,914 placements, 4 wrong (0.008%).
  - O2, on the reads within 10 AS of a rival: 7,000 reads. It placed 349, all right, and abstained on 6,651. v14 decided the
    same 7,000 identically.
  - 216 reads were `as_unexplained` (left to O3), 5 `as_elsewhere` and 32 not hosted.
- **The cost runs' own correctness bars** (the mainstay set, development):
  - A119b passes.
  - public fails one bar: INTACT B2, 25 errors against 24.

**What it means:**
- **The split holds where O2 acts.** On the AS-tied reads, O2 is never wrong, as v14 was not.
- **The split's errors are minimap2's.** It is wrong on 4 of its 49,914 resolved placements. Those reads' best AS lies on
  the wrong copy by 10 or more, here an unspliced alignment to a processed pseudogene copy.
- **Off-family ties make the difference** between v16 and v14: 35 fewer errors there, with v16 abstaining on 348 fewer
  correct placements.
- **Not adopted as registered.** The B2 bar is strict, as the user chose, and one read failed it in three worlds.


---

## PREREG_o3_graph_round22_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_graph_round22_2026-10-19.md`)*

## Round 22: O2 only on ambiguous reads, with sequence ties, on fresh multigene families (2026-10-19)

This pre-registration was written before any world of round 22's families was built; no read of these families has been
simulated. Code: branch `o3-fast`, frozen at the commit that adds this file.

### Why

**Round 21 failed by one read** (register 1331): 5 of 6 cells passed, with C1 passing.
- **The read:** NCOR1P4. Its own copy aligns all 340 bases across 3 introns (AS 262, NM 0). The higher-scoring copy leaves
  38 bases soft-clipped and has 1 intron (AS 278, NM 0).
- **Why it was misplaced:** minimap2's AS charges every intron and nothing for a soft clip, so the margin of 16 made the read
  minimap2's, and it went to the copy that explains less of the read.

**The user (2026-10-19), on the proposed fix:** "ok, lets do that". A read is ambiguous, and so O2's, when its AS is within
10 of the best (Eichler's rule), or when an alignment on another locus explains it with no more edits than the best-AS one,
counting every base it leaves unaligned (`veto_distance`).

**Development (all spent):**
- **Round 21's 36 worlds:** 6 of 6 cells now pass.
- **Round 19's:** 4 of 6, unchanged. The failures are the PABPC1L2B reads, which fit the other copy better by sequence
  (NM 1-2 against 7-8): the split's accepted cost, under the strict B2 the user chose.
- **HG002's untied reads with a donor copy:** 44,086 right and 594 wrong, against v14's 44,021 and 599.
- **HG002's tied reads:** 50 and 7, as v14.

### The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v17** = v16 (round 21) + `--seq-tie`. That is v14 + `--o2-scope tied --as-slack 9 --seq-tie`
(`o2graph.py run W intega`).
- **A read's loci:** its alignments in `E0.bam` within 9 AS of its best, plus any alignment whose `veto_distance` is at
  most the best-AS alignment's.
- **The decision rules are round 21's:**
  - one seeded copy: minimap2's placement, unless its alignment does not explain the read;
  - a tie: v14's call when it is one of the tied copies;
  - a locus off the family's copies: abstain;
  - the silent copy alone: abstain.

**Integrated flags (`$V17`):** round 21's `$V16` plus `--seq-tie`:
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
--o2-scope tied --as-slack 9 --seq-tie
```

**The simulation:** as rounds 19-21 (profiles v4).

### Families (held out; drawn by the fixed rule)

**The pool:** the adopted multigene rule with all components (`r19seq/families19cal_full2_all.json`, 439 families).
- **Spent:** every family touching a gene of an earlier simulated set (27 files: `ls families*.json | grep -v _all.json`,
  round 21's included), or of the HG002 worlds.
- That leaves 98 unspent families with 303 copies. File: `r22seq/families22seq_all.json`.

**The draw:** the rule of rounds 17-21 with seed 22.
- **26 families, 90 copies:** 38 coding and 52 non-coding (34 transcribed processed, 16 transcribed unprocessed,
  1 unprocessed, 1 IG_V).
- **Copies per family:** {2: 19, 3: 3, 6: 2, 8: 1, 23: 1}.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families22.json` (set `r22`, sha1 4aadf1d72ccc).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19ABCC6P2 | 2 | - | ABCC6P2 (tr. unproc.), ABCC6P1 (tr. unproc.) | chr16 |
| S19AC1143242 | 2 | - | AC114324.2 (tr. unproc.), AC114324.2 (tr. unproc.) | chr5 |
| S19AC1364281 | 2 | AC136428.1 | IGHV3-11 (IG_V_gene) | chr14, chr16 |
| S19AC2337026 | 23 | - | AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.) | chr17 |
| S19ANKRD18B | 2 | ANKRD18B | ANKRD18CP (unproc.) | chr9 |
| S19ANKRD30BL | 6 | ANKRD30BL | FO681491.1 (tr. unproc.), CR383658.1 (tr. unproc.), AP005901.1 (tr. unproc.), ANKRD30BP2 (tr. unproc.), CR383658.1 (tr. unproc.) | chr2, chr13, chr14, chr18, chr21, chr22 |
| S19CDK2AP2P1 | 8 | - | CDK2AP2P1 (tr. proc.), PTGER4P2 (tr. proc.), CDK2AP2P1 (tr. proc.), CDK2AP2P1 (tr. proc.), PTGER4P2 (tr. proc.), PTGER4P2 (tr. proc.), PTGER4P2 (tr. proc.), PTGER4P2 (tr. proc.) | chr9, chr13, chr14, chr21 |
| S19CYP11B1 | 2 | CYP11B1, CYP11B2 | - | chr8 |
| S19DEFB105A | 3 | DEFB105A, DEFB105B | - | chr8 |
| S19FAM205A | 2 | FAM205A | FAM205BP (tr. unproc.) | chr9 |
| S19FOLH1 | 2 | FOLH1 | FOLH1B (tr. unproc.) | chr11 |
| S19GPAA1P1 | 2 | - | GPAA1P1 (tr. proc.), GPAA1P2 (tr. proc.) | chr2 |
| S19KRT17P1 | 3 | - | KRT17P1 (tr. unproc.), KRT17P2 (tr. unproc.), KRT17P6 (tr. unproc.) | chr17 |
| S19LDHAL6B | 2 | LDHAL6B | LDHAL6CP (tr. proc.) | chr12, chr15 |
| S19LILRA2 | 2 | LILRA2, LILRA1 | - | chr19 |
| S19MT1G | 2 | MT1G, MT1H | - | chr16 |
| S19PCDH11X | 2 | PCDH11X, PCDH11Y | - | chrX, chrY |
| S19PDXDC1 | 2 | PDXDC1 | PDXDC2P (tr. unproc.) | chr16 |
| S19PRIM2 | 2 | PRIM2 | AC244258.1 (tr. unproc.) | chr6 |
| S19SPIN2B | 2 | SPIN2B, SPIN2A | - | chrX |
| S19TYW1 | 2 | TYW1, TYW1B | - | chr7 |
| S19ZNF208 | 3 | ZNF208, ZNF676, ZNF728 | - | chr19 |
| S19ZNF506 | 6 | ZNF506, ZNF93, ZNF486, ZNF737, ZNF626, ZNF66 | - | chr19 |
| S19ZNF611 | 2 | ZNF611, ZNF600 | - | chr19 |
| S19ZNF799 | 2 | ZNF799, ZNF443 | - | chr19 |
| S19ZNF98 | 2 | ZNF98, ZNF492 | - | chr19 |

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 1201, 1202 and 1203, with no seed selection. That makes 6 cells and 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families22.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout22
for S in 1201 1202 1203; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r22 $S --local $O/local $V17 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

### Bars

**Correctness:** as round 21 (fastset's own e2e step).
- **Table and scoring:** the scoped table `o2.intega.m2.tsv`, per-read detectability less each world's error floor.
- **Set apart:** reads tied with a locus on no copy of the world (`--off-family-ties 9 --seq-tie`) are the class
  `off_family_tie`, left out of both systems' sums.
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**Verdict.** A cell PASSES when all six of its worlds pass both bars, and round 22 PASSES only when all six cells pass (72 bar
checks).

**Cost (C1),** as round 21, on the mainstay set (seed 306, both profiles), into `fast/r22cost/v17/`.
- **Bar:** v17's integrated CPU seconds (O1 + O2 + O2_genome + O2_scope + O3) at most the current method's, in each profile.

**Reported, no bar:**
- v14 on the same worlds and scoring.
- The off-family ties' fates.
- The split of placements between minimap2 (`as_resolved`) and O2: O2's domain size, its placements and errors.
- Wall times and peak RSS.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

### What a PASS would and would not claim

**A correctness PASS** would claim that O2 confined to ambiguous reads, with ambiguity defined by AS within 10 or by
equal-or-better sequence on another locus, does at least as well as the current method in every world type. It would hold
on 26 held-out multigene families, most of them with processed pseudogene copies.

**It would not claim:**
- **Anything about reads tied with a locus off every copy.** These are scored apart.
- **Accuracy where a donor's copy carries its paralog's alleles.** Its reads fit the paralog better by sequence; register
  1330 and round 19's development.
- **Real-read accuracy beyond the HG002 development numbers above.**

### Amendment 1 (2026-10-19, 11:25, before any crashed run was rerun)

**The crash.** `/mnt/linuxdisk` filled at 11:17, while cell public_s1201 was writing its masked genome.
- The other four cells and both cost runs then failed at their first write. Their output files are empty.
- **Cell A119b_s1201 had already completed: PASS, all six worlds.**

**The fix.** The user approved deleting the masked genomes of the finished rounds (regenerable from each cell's seed) and
the day's HG002 temporary files: 32 files, 86 GB. The six read-only round directories, rounds 9-14, were left alone.

**The rerun.** The crashed runs are rerun once, unchanged, as the rules require: public_s1201 (its partial `masked/` and
`set.log` cleared), A119b_s1202, public_s1202, A119b_s1203, public_s1203, and both cost runs.

### Results (2026-10-19; `/mnt/linuxdisk/tmp/o3_graph/fast/heldout22/`)

**Correctness: PASS.** All six cells pass and all 72 bar checks hold: A119b_s1201 in the first run, the other five in the
Amendment-1 rerun, with no crash.

| system | judged correct | judged errors | off-family ties: correct / errors |
|---|---|---|---|
| integrated v17 | 78,421 | 3,255 | 3,355 / 94 |
| current | 53,763 | 12,600 | 1,652 / 380 |

**C1: incomplete.**
- A119b passes: 953 integrated CPU s against the current method's 1,086.
- The public cost run was stopped by the user at 12:35 (2026-10-19), before it finished, to refocus on the families where
  O2 can be shown to work. The round's cost bar is therefore UNDECIDED.
- Rounds 19-21 passed C1, and the scope step adds about 5 s.

**What it means.** With ambiguity defined by Eichler's within-10 rule plus sequence ties, O2 confined to the ambiguous reads
does at least as well as the current method in every world type, on 26 held-out families that are mostly processed
pseudogenes. Adoption as the default needs the cost bar's public profile.


---

## PREREG_o3_npip_candidates_realign_2026-10-21.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_npip_candidates_realign_2026-10-21.md`)*

## Pre-registration: do the maternal-only NPIP locus's reads find their copy once the O3 candidates join the reference? (2026-10-21)

**User (2026-10-21):** "ok lets do that". The proposal: add the two candidates from register 1363 to the reference as extra loci,
realign, and see whether C's 37 reads move onto its candidate. This is the read-recovery step of the 10-01 chain (Amendments 5-9
of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md`, arms R / R+I+L / C), applied to register 1363's output.

**Question.** IsoCon left most of the locus's reads unclustered (register 1363: C 7 of 37). Does adding the candidates' sequences to
the reference give the reads a home through the aligner, without pulling in other genes' reads?

### Fixed inputs (from register 1363, frozen in `/mnt/linuxdisk/tmp/o3_npip/audit_frozen_4c64ec81/`)

- **Candidates added:** the two counted candidates (>= 2 contigs).
  - **B** = `iso_NPIP_14` + `iso_NPIP_15` (his maternal chr18 6.744-6.756 Mb);
  - **C** = `iso_NPIP_16` + `iso_NPIP_17` (6.767-6.804 Mb).
  - The single-contig `iso_NPIP_19` is under the floor and is not added.
- **Reads:** the 2,955 net reads.
- **Units** (the audit's definition, fixed here before any realignment). A target read (register 1363's R_T, 94 reads) is labelled
  by the aligned bases (>= 50) of its best alignment in T on his maternal chr18:
  - **A** = 6,720,000-6,744,300; **B** = 6,744,300-6,756,100; **C** = 6,766,000-6,822,000;
  - so a read is C, B, A, or "B into A".
  - The audit counted C 37, B 21, B into A 21, A 11.

### Arms

Both arms use the library's command, `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (minimap2 2.30).
- **R0:** the net reads realigned to GGO alone (`GGO.splice.mmi`).
- **C:** the same reads aligned to GGO plus the four contigs (`minimap2 -x splice -d`, as Amendment 9 built its index).

A read's call is the locus of its primary: a candidate (B or C, either of its contigs) or a reference position.

### Bars

- **H1 (C's reads find their copy):** in arm C, at least half of unit C's reads have their primary on C's candidate.
- **H2 (B's reads):** at least half of unit B's reads have their primary on B's candidate.
- **H3 (no false moves):** at most 1% of the non-target net reads (2,955 − 94) have their primary on either candidate in arm C.

**PASS** = H1, H2 and H3; **FAIL** otherwise, naming the bars that fail.

**Reported, not bars:**
- the "B into A" and A reads;
- register 1360's 29 reads;
- the MAPQ of reads that moved;
- the target reads AS-tied (within 10) among NPIP copies in R0 vs C, i.e. the O2 domain this O3 step would shrink.

### What a PASS would claim, and what it would not

**Claim.** In his own data, a copy the reference lacks, once named by the reference-free chain, takes back its reads through ordinary
alignment, at a cost of at most 1% of other reads.

**Not claimed:**
- that the chain finds such copies in general (register 1363 is one locus);
- unit A (no candidate);
- other samples.

### Rules

- Run once; bugs → amendment + rerun.
- The index build and the alignments run in the foreground, one at a time, each under 10 minutes.
- The metric-auditor checks the verdict before it is reported.
- Gorilla only. Push only when asked.
- Driver: new `realign` and `score-realign` steps in `bench/rna_allele/npip_o3.py`, with tests.

### Results (driver 5d04d56e; run once; `/mnt/linuxdisk/tmp/o3_npip/realign/`): **PASS**

- **Index:** GGO + the 4 contigs, built in 3 min 26 s. Both arms aligned all 2,955 reads (R0 2 min 24 s, C 2 min 18 s).
- **Units as registered** (no strand): C 39, B 23, A 11, B into A 21. The audit's 37 / 21 were the main-strand reads; 2 C and 2 B
  reads lie on the opposite strand.

| bar | arm C | bar value | holds |
|---|---|---|---|
| H1: unit C's reads on C's candidate | **33 / 39** (85%; main strand 33 / 37) | >= 20 | yes |
| H2: unit B's reads on B's candidate | **21 / 23** (91%; main strand 21 / 21) | >= 12 | yes |
| H3: other reads moved onto a candidate | **0 / 2,861** | <= 28 | yes |

**Reported:**
- **Where each unit's reads go in arm C:**
  - C: 33 on C's candidate, 6 on the reference;
  - B: 21 on B's candidate, 2 on the reference;
  - A: 6 on C's candidate, 5 on the reference;
  - B into A: 2 on B's candidate, 19 on the reference.
- **The 4 opposite-strand reads** all stay on the reference.
- **Register 1360's 29 reads** (24 in unit C, 5 in A): 24 on C's candidate (21 C reads + 3 A reads), 5 on the reference (3 C
  reads + 2 A reads).
- **MAPQ of target reads on a candidate:**
  - C: 26 at >= 30, 10 at 1-29, 3 at 0;
  - B: 23 at 1-29 (its two contigs are near-identical);
  - homed reads split between C's contigs (iso_NPIP_16 18, iso_NPIP_17 15).
- **The O2 domain:** target reads AS-tied within 10 among loci fall from **34 (R0) to 10 (C)**.

**Reading.** Once the chain names the missing copies and they join the reference, ordinary alignment gives their reads a home:
89-100% of the main-strand reads of the two named units, with no other read pulled in. Reads of unit A, which got no candidate, stay
on the reference, though 6 of its 11 land on C's candidate. Most of the tied reads that O2 would have faced at this locus stop
being ties.

Audit: below (PASS holds with caveats; the reading above is corrected).

### Audit (metric-auditor, 2026-10-21): PASS holds with caveats

Every number reproduces:
- **Units:** C 39, B 23, A 11, B into A 21.
- **Bars:** H1 33 / 39 (main strand 33 / 37); H2 21 / 23 (21 / 21); H3 0 / 2,861.
- **Reported figures:** the per-unit fates, register 1360's split (21 + 3 / 3 + 2), the MAPQ counts and the ties (34 → 10).
- **Order holds:** the design came before `units.tsv`, the index, R0, C and the results.
- **Strand:** the strand-free labels were applied as defined, and the main-strand reading passes as well. The 37 / 21 / 21 / 11
  quoted in the design do not follow from its own definition.

**Caveats:**
- **H1 is robust.** On reads IsoCon did not use, 26 / 32 (81%). Under any tie policy it is 0.83-0.89.
- **H2 is mostly self-fulfilling.**
  - 19 of B's 21 homes were IsoCon members, the very reads the candidate was built from. On non-members H2 is 2 / 4, exactly at the
    bar.
  - The "B into A" reads start where B's do. Counted as B's reads, B's candidate holds 23 / 44, and 4 / 25 non-members.
- **B is not an NPIP copy.** No B read has an aligned block on an NPIP copy. B's reads entered the net only because introns of
  >= 100 kb span NPIP.10 (the net collects reads by span overlap).
- **H3 could not fail.** Only 16 of the 2,861 non-target reads have any record on a candidate, below the bar of 28; the closest
  scores 735 AS below its primary.
  - Reads of B's relatives were never realigned: RRN3 has 1,019 reads in the library (31 in the net), the RRN3-like LOC101142371 65,
    a locus at 29.42 Mb 4,764.
  - "No other read pulled in" therefore holds for NPIP-net reads only.
- **Unit A is C's inverted duplicate.** 10 of 11 A reads fit C equally well, 8 of them with exactly the same `ms`. The 6 "A reads
  on C's candidate" are exact ties, so their origin is undecidable; the same holds for register 1360's 5 A reads.
- **The arms differ beyond the contigs.** 26 of the 2,845 non-target reads with no candidate record changed primary between R0 and
  C. Comparisons between the two arms carry that noise; the bars do not.
- **Reported definitions:**
  - The ties (34 → 10) count any locus, 100-kb bins included. The design's wording (ties among NPIP copies) gives 33 → 8.
  - B's MAPQ 1-29 comes from competition between its own two contigs: the other contig is the runner-up for all 23 reads, and in
    R0 the same reads had MAPQ >= 30.
- **Not held out.** 10-01's Amendment 5 already moved this locus's reads onto IsoCon transcripts. This is a replication at one
  locus, with no oracle arm.

**What stands:** PASS as registered, with the weight on H1. Adding the reference-free candidate for register 1360's copy gives 33
of its 39 reads a home (26 of 32 among reads the candidate was not built from), and moves no other NPIP read. Whether reads of B's
relatives in the RRN3 family would be pulled in is untested.

### Amendment 1 (2026-10-21, before any of its alignments): the candidates' relatives (user: "Yes please proceed")

**Why.** The audit found H3 inert: only 16 of the 2,861 net reads could reach a candidate at all. Reads of the candidates' own
sequence relatives were never realigned, because they lie outside the NPIP net.

**Read sets** (fixed by `realign/rel.names` and `rel.sets.tsv` before any alignment): library reads with any record at these loci
and not in the net.
- **RRN3 family** (B's relatives, the audit's three loci on NC_073242.2): 5,863 reads in all.
  - RRN3, 29,559,520-29,593,927;
  - LOC101142371, 31,545,000-31,600,000;
  - the 29.42-Mb locus, 29,415,000-29,425,000.
- **C's relatives:** the aligned blocks of C's contigs' GGO hits (`outputs.pri.paf`) outside the 25 NPIP copies, 10 blocks:
  2,476 reads.
- **Union:** 8,339 reads, each from its primary record in sequencing orientation.

**Arms.** R0 (GGO) and C (GGO + the 4 contigs), with the registered command, 2 chunks per arm.

**Bar H3′:** in arm C, at most 1% of the 8,339 reads (83) have their primary on a candidate.
- **Reported, not a bar:** per set; the reads reachable (any record on a candidate); the near-tied (a candidate alignment within 10
  AS of the read's best).

**Consequence.** H1-H3 and the registered PASS stay as they are.
- If H3′ fails, "no other read pulled in" is withdrawn.
- If it holds, that statement extends to the candidates' sequence relatives.

### Results under Amendment 1 (driver 05fdaed4; run once): **H3′ holds**

- **Alignment:** all 8,339 relatives' reads aligned in both arms, 2 chunks each, about 2 min per chunk.

| set | reads | arm R0: on a candidate | arm C: primary on a candidate (moved) | reachable (any candidate record) | near-tied (within 10 AS) |
|---|---|---|---|---|---|
| RRN3 family | 5,863 | – | **1** | 86 | 1 |
| C's relatives | 2,476 | – | **0** | 1 | 0 |
| all | 8,339 | – | **1** (bar <= 83) | 87 | 1 |

- **The test could fail this time:** 87 reads reach a candidate (against 16 in H3's net), and 1 of them moves.
- **The moved read** (SRR27438212.2412175, RRN3-family set):
  - in R0 its primary is a poor fit: NC_073242.2:29,462,744, AS 370, its secondaries at 32.19 Mb reaching 497;
  - in C it fits B's candidate at AS 725 (MAPQ 19; the runner-up is B's other contig, at 712).
  - So it fits B far better than anything in the reference: possibly a B read the span-based net missed. This is not verified
    against his assemblies.

**Reading.** Adding the two candidates pulls in no read of their sequence relatives outside the NPIP net (at most 1 of 8,339, and
that one fits B better than any reference locus). The no-false-move statement now extends to the candidates' relatives.

Audit: below (H3′ holds with caveats; the reading above is corrected).

### Audit of Amendment 1 (metric-auditor, 2026-10-21): H3′ holds with caveats

Every count reproduces (moved 1 / 1 / 0, reachable 87 / 86 / 1, near-tied 1 / 1 / 0 for all / RRN3 / C's), with 8,339 primaries per
arm, one per read. The read set rebuilds exactly and was fixed before the alignments.

**Corrections:**
- **The mover is a target read, not a false move.** 97.7% of its bases are exact 25-mers of his maternal chr18, B
  (6,745,103-6,756,064) running into A: a "B into A" read that the span-based net missed. False moves: 0.
- **Most of the set are passengers.** Only 951 of the 8,339 reads have an aligned base on a candidate's GGO block (947 in the RRN3 set,
  4 in C's). The 29.42-Mb window's 4,742 are exon reads of another gene (LOC101140878), and 2,472 of C's 2,476 only span a block with
  an intron. C's relatives' blocks are 4 distinct loci, not 10 (the set is unaffected).
- **"Could fail" holds only in arithmetic.** Every non-mover sits >= 435 AS below its best. The one near-tied read is the mover
  itself.
- **Missed target reads:** 5 more reads in the set match the maternal locus at >= 93% and stayed on the reference. They are likely
  target reads the net missed, not verified against his assemblies.

**What stands:** adding the two candidates moves none of the 951 reads aligned at B's RRN3-family relatives' blocks (the one mover
is a target read). The no-false-move statement extends to B's RRN3 relatives. C's NPIP-like relatives outside the net are untested
(4 reads).


---

## PREREG_o3_npip_maternal_copy_2026-10-21.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_npip_maternal_copy_2026-10-21.md`)*

## Pre-registration: does the reference-free chain name KB3781's maternal-only NPIP locus as exactly one new copy, blind to his assemblies? (2026-10-21)

**User (2026-10-21):** "ok lets do that". The proposal: cluster his NPIP reads, link alleles at 0.96%, and check that exactly one new
candidate appears and that it holds the 29 reads of register 1360. His assemblies serve only as the answer key.

**Question (O3: detect and flag).** Run on KB3781's own fibroblast reads of the NPIP family, with nothing masked, does the frozen
reference-free chain:
- name exactly one candidate new copy, namely the maternal-only locus whose reads minimap2 forces onto reference NPIP copies;
- hold that locus's reads;
- raise no false flag?

### The chain (frozen)

The chain is that of Amendments 7-9 of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md`, with code in
`bench/rna_allele/link_test.py`, `merge_test.py` and `control_test.py`. It never sees his haplotype assemblies.

1. **Net.** Every read of his fibroblast library (`fibroblasts/GCA_029281585.2_flnc_mm.bam`, aligned to GGO with `-N 50 -p 0.1`)
   with any alignment record on one of the 25 NPIP copies (catalog family `NPIP`). Each read's sequence comes from its primary record,
   in sequencing orientation. Two deviations from Amendments 7-9:
   - **No 1,000-read cap.** The cap only bounded runtime: IsoCon ran the 2,951-read NPIP net in 93 s on 10-01, and a cap would
     subsample the target's reads.
   - **No unmapped reads.** His 959 unmapped reads have a median length of 69 bp, and 10-01 found no family evidence in them.
2. **IsoCon** 0.3.3: `IsoCon pipeline --nr_cores 4`, defaults.
3. **Flag.** An output is "not in the reference" when its best spliced alignment to GGO (`minimap2 -c -x splice:hq -uf -N 20`, the
   GGO splice index) has identity × coverage < 0.999.
4. **Link.** A flagged output whose whole-length divergence to that best locus is d <= δ = 0.00958 is an allele of it; otherwise it
   is a new-copy contig. δ is the 99th percentile of exonic divergence between his two haplotypes over 28,541 single-copy genes,
   frozen on 10-01.
5. **Merge.** New-copy contigs are aligned all-vs-all (`minimap2 -c -x asm20 --cs -N 200 -p 0.1`). Two join when their best
   alignment covers >= 50% of the shorter contig and de <= δ. The components are the candidate copies.
6. **Support floor.** A candidate counts only with >= 2 contigs. Amendment 9's no-deletion control proposed this floor, after the
   fact, for the next tests; it is fixed here before any run. One-contig candidates are listed, not counted.
7. **A candidate's reads:** the reads IsoCon used for its contigs (`cluster_info.tsv`).

### Answer key (his assemblies; evaluation only)

- **Candidate class** (Amendment 9's `classify`): the best contig's spliced alignment to his maternal and paternal assemblies
  (`-c -x splice:hq -uf -N 20`, the haplotype splice indexes).
  - identity × coverage < 0.999: **unmatched** (false);
  - on the haplotype the reference is identical to: **pri** (false);
  - inside the lift of an NPIP copy (the frozen 10-01 asm5 lift, `truth_lift.py`, copy spans): **allele** (false);
  - otherwise: **haplotype-only** (true).
- **Target T:** the maternal-only locus of register 1360, his maternal chr18 `CM054600.2:6,720,000-6,820,000`.
- **Target reads R_T:** net reads whose best alignment in his genome lies in T and scores >= 5 `ms` above their best alignment
  outside T. Alignments are `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1` to each haplotype assembly, compared by `ms`;
  alignments with an intron of >= 100 kb are dropped. The 29 tied reads of register 1360 are reported among them.

### Bars

- **E1 (found):** at least one counted candidate is haplotype-only, with its hit inside T.
- **E2 (one copy, one candidate):** exactly one such candidate.
- **E3 (holds its reads):** that candidate holds at least half of R_T.
- **E4 (pure):** at least 90% of its reads are in R_T.
- **E5 (no false flag):** no other counted candidate is false (unmatched, pri or allele). Other haplotype-only candidates are
  reported as additional flags, not failures; one example would be 10-01's maternal-only copy at chr18 mat 17.40 Mb.

**PASS** = E1-E5. **FAIL** otherwise, and the report names the bars that fail.

### Seen before: this is a positive control of the chain, not a blind discovery

- Register 1360 found the locus from his assemblies: 29 tied NPIP reads with 0-3 substitutions there and 18-167 at their best
  reference placement.
- On 10-01, IsoCon on the same library (with NPIPA2 masked) produced 4 outputs from this sequence (chr18 mat 6.72-6.76 Mb, 79-93% to
  the reference), and Amendment 5 moved 31 reads onto two of them (iso_31 / iso_32).
- **Never run before:**
  - the full chain (flag, link, merge, support floor) on NPIP with nothing masked;
  - its read accounting against R_T;
  - its false flags within NPIP.

### Rules

- Run once; bugs → amendment + rerun.
- Heavy commands in the foreground, one at a time, each under 10 minutes.
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Gorilla only. Push only when asked.
- Code: a new driver, `bench/rna_allele/npip_o3_test.py`, calls the frozen pieces (the flag and link of `control_test.contigs`,
  `merge_test.pairs` / `components`, `truth_lift`), with tests for its own new parts (the net, R_T and the bars).

### Results (driver a7a81ece, `bench/rna_allele/npip_o3.py`; run once; `/mnt/linuxdisk/tmp/o3_npip/`): **FAIL** (E3, E4, E5)

- **Chain:**
  - net 2,955 reads (10-01: 2,951);
  - IsoCon: 109 outputs, 1,600 reads used (54%);
  - 21 flagged not in the reference; 16 linked back as alleles; 5 new-copy contigs;
  - **3 candidate copies**, 2 counted (>= 2 contigs).
  - The 25 copies all lift (T2d 16, T2i 6, T? 3).
- **Target reads R_T:** 94 net reads fit his maternal chr18 6.72-6.82 Mb best by >= 5 ms. All 29 of register 1360's are among them.

| candidate | contigs | best fit in his genome | class | reads (in R_T) | best reference locus |
|---|---|---|---|---|---|
| iso_NPIP_16 + 17 | 2 | mat chr18 6,766,688-6,803,775, 0.9991 | haplotype-only, **in T** | 8 (7) | NPIP.1 |
| iso_NPIP_14 + 15 | 2 | mat chr18 6,744,425-6,756,028, 0.9986 | unmatched (< 0.999), **in T** | 19 (19) | NPIP.10 |
| iso_NPIP_19 | 1 (not counted) | mat chr18 100.86 Mb, 0.8718 | unmatched | 5 (1) | elsewhere |

- **E1 holds, E2 holds:** exactly one counted haplotype-only candidate lies in T.
- **E3 fails:** it holds 7 of the 94 target reads (bar 47).
- **E4 fails:** 7 of its 8 reads are target reads (0.875; bar 0.90).
- **E5 fails:** the other counted candidate is "unmatched", so it counts as false by the frozen class. But it lies in T as well, 0.0004
  below the 0.999 consensus bar.
- **Where the 94 target reads went:** IsoCon left **66 unassigned (27 of register 1360's 29)**; 26 are in the two target candidates
  (pure: 26 of their 27 reads); 1 is in an output linked back as an NPIP.5 allele; 1 is in the uncounted candidate.

**Reading.** The chain raises no flag outside the maternal-only locus, apart from one single-contig candidate under the floor. It
does find that locus, but splits it in two and recovers only 28% of its reads. The loss is upstream of the flag-link-merge steps:
IsoCon does not converge on most of the locus's reads, as on 10-01 ("real reads too truncated and sparse").

Audit: below (FAIL holds with caveats; the reading above is corrected).

### Audit (metric-auditor, 2026-10-21): FAIL holds with caveats; the reading above is corrected

Every number reproduces:
- **Chain:** net 2,955; 109 outputs; 1,600 reads assigned; 21 flagged / 16 linked / 5 new-copy contigs; 3 candidates.
- **Bars:** E1-E5 = holds, holds, fails, fails, fails.
- **R_T:** 94 reads. The same set comes out at ms margins of 0, 1 and 10; without the intron filter it is 93. Every target read fits T
  best by >= 103 ms.
- **Fates:** 66 / 19 / 7 / 1 / 1. Register 1360's reads: 27 of 29 unassigned.
- **Code path:** only the two declared deviations. The code is byte-identical to the results commit, and the chain is unchanged since
  10-01 apart from docstrings.

**Corrections:**
- **T holds three transcription units, not one copy split in two.**
  - C (6.767-6.804 Mb) is register 1360's copy. It holds 24 of the 29 reads; the other 5 are in A (6.72-6.74 Mb), none in B.
  - B (6.744-6.756 Mb) is a separate unit on the same strand, with its own start site: 21 reads begin at 6,756,009-6,756,064. Its 5′
    exons match RRN3's first three exons (GGO chr18 29.59 Mb). No read joins B and C, and their contigs share no sequence.
  - **Each unit that got a candidate got exactly one:** C → iso_NPIP_16 + 17 (haplotype-only), B → iso_NPIP_14 + 15. A got none.
- **E5's "false" flag is consensus noise on a real unit.** iso_14 and iso_15 each differ from his assembly by ±1 A in poly(A) runs of
  16-17 nt, and the reads covering those runs disagree too (26 / 36 and 22 / 36 carry indels).
- **E1 and E2 hold by one base.**
  - iso_16 scores 1055 / 1056 = 0.99905; iso_17 alone would be unmatched at 0.99898.
  - Had iso_15 matched one base better, E2 would fail instead of E5.
  - The bars could not credit a second true unit in T: a design flaw, since the target window was assumed to hold one copy.
- **E3's bar was out of reach even for a perfect candidate.** No single unit has 47 reads (C 37; B 21, plus 21 that run on into A).
  Per unit:
  - C: 7 / 37 held, 30 unassigned;
  - B: 19 / 21 held;
  - B running into A: 0 / 21;
  - A: 0 / 11.
- **E4 is genuine, but on n = 8.** SRR27438212.5855649 has no record in T and fits pat chr18 22,085,194 at NM 0.
- **Table and text:**
  - "NPIP.10" for iso_14 + 15 came from a 172-kb GGO hit with a 168,474-bp intron; none of its exons touches NPIP.10.
  - "No flag outside T" should read "no new-copy candidate outside T": 16 flagged outputs were linked back.
- **History:** IsoCon converges where a unit's 5′ ends agree (B: 19 / 21) and fails where they spread (C spans 44 kb: 7 / 37). That
  matches 10-01 for C only.

**What stands:** FAIL as registered.
- The chain raised no new-copy candidate outside the maternal-only locus.
- It named the locus's two transcribed units with one candidate each, one of them misclassed by a poly(A) homopolymer error.
- It recovered only 7 of the 37 reads of register 1360's copy.

Detection works here; read recovery depends on IsoCon converging, which fails on transcripts whose 5′ ends spread.


---

## O3_FINGERPRINT_SURVEY_2026-10-22.md

*(original: `docs/archive/2026-10/merged/O3_FINGERPRINT_SURVEY_2026-10-22.md`)*

## O3 as an RNA-only "look harder" screen: synthesis of the vault survey (2026-10-22)

**Status:** synthesis of four survey files (`docs/o3survey_2026-10-22/`, angles A–D) plus two quick checks run on the frozen O2-VC
outputs. It is design and exploration only. No screening score was computed on held-out truth, so HG00733 and HG002 stay usable for
evaluation.

### The framing

- **The advisor's view (Canzar).** "Cluster by PSV, BAM tags and find a certain print of a missing copy"; it "does not need to be
  definitive (that is what DNA WGS is for) but should suggest when to look harder".
- **The user's rule (2026-10-22):** RNA reads against one reference only. No DNA reads, no per-sample assembly.

So O3 is a **screening rank**, judged by AUC or by sensitivity at a fixed false-alarm rate. It is not a hard abstention gate for O2.

### Two statements that bound what is possible

**1. An extra copy can be certified from RNA reads** (survey C, statement 1; survey A, fingerprint 1).

- **The test.** Take the read patterns over a family's columns, where each pattern is seen in ≥ 2 reads and two patterns conflict at
  ≥ 2 columns. If the patterns that mutually conflict cannot all fit in the 2 allele slots per copy of the copies the reads tie to,
  the donor carries more sequence than the reference: an extra copy.
- **Properties.** The test is exact, checked by bipartite matching, and uses no thresholds. The "≥ 2 independent observations"
  rule handles noise. Comparing the haplotype count with 2 × the copy count is its weak special case.

**2. An absent copy cannot be certified from RNA alone** (survey C Theorem A; survey A's verdict; survey D's register evidence).

- **The argument.** Suppose the donor carries copy X as h1/h2 and lacks Y. Another donor carries X as h1/h1 and Y as h2/h2. Both
  express the same read sequences inside the shared block, so every RNA test behaves the same on both.
- **What follows.**
  - The PRAMEF7.1 case can only ever be a "consider" flag.
  - H5's lost right placements are forced, not a tuning failure.
- **The single exception.** A read that reaches the copy's own unique sequence (a flank, a unique exon, a readthrough) certifies
  that the copy is **present**. That is survey A's "anchored ownership".

### The screen: four tiers, threshold-free (ranked by tier, then by consequence)

| tier | fingerprint | what it means | source | status |
|---|---|---|---|---|
| 1 certified extra copy | haplotype excess: mutually conflicting patterns exceed the tied copies' allele slots | the donor has a copy CHM13 lacks | C1, A1, B1 | needs each copy's **teacher** alleles. Re-call with longcallR's haplotagged BAM (we ran `--no-bam-output`) |
| 2 consider: novel recurrent haplotype | a recurrent haplotype no copy's teachers show; generalises the current linked foreign-allele flag to whole haplotypes | extra copy, mosaic/fusion, or an untaught allele | C3, D3 | the current flag already fires on HG002's NPIPB2 (donor copy with a 57-bp insertion), SPDYE1 (4 flags; 28 real donor-insertion reads) and the CNV families NCF1, NOMO1, GOLGA8A, SMN2. It does not fire in HG00733's 46 absent-copy families (0/46) and false-alarms in 5 of 128 |
| 2b consider: unexplained recurrent indel | an indel no tied copy explains; a copy just outside the 10-AS tie has fewer edits | a structurally different donor copy (NPIPB2.2) | B2 | screen; also fires on simulated KRTAP21 |
| 3 consider absent | collapse-compatible pair (two copies' reads fit in two haplotypes) AND no read anchors either copy to its own unique sequence | a reference copy the donor may lack (PRAMEF7.1) | C4, A2, A3, H5 | necessary, never sufficient (statement 2). Its false-alarm set (present homozygous pairs) is computable on the truth side before running |

**Dead, do not re-propose** (survey D and the register):
- MEC used as a haplotype count;
- colouring the raw conflict graph;
- the four-gamete test (it never certifies when ≥ 2 copies are tied);
- depth or expression as a copy-number signal;
- soft clips;
- raw NM / `de` (HiFi indel error dominates);
- reference-PSV-only O3.

**Prior evidence for tier 2's lineage** (survey D):
- the linked-cluster to IsoCon chain recovered 425/426 in simulation with 0/152 false;
- real data were weaker (r1292–1296: not claimed; r1304: alleles falsely flagged in 4 of 8 intact families);
- IsoCon detected both of NPIP's maternal-only units with nothing outside (r1363).

### Quick checks run tonight (descriptive, post hoc)

- **The current flag as a family-level screen:**
  - HG00733: 0/46 absent-copy families flagged, 5/128 others flagged. That is expected, since a reference copy the donor lacks
    creates no foreign alleles, and HG00733's simulation contains no donor-only copies.
  - HG002: it fires in the families with real structural differences listed in the table.
- **A family-level haplotype-excess check on tied reads only** (survey A): 0 of 131 families in HG00733 and 0 of 50 in HG002, so
  specific but blind. That includes SPDYE1 and RASA4B, which do contain real donor insertions. The certificate needs the **teacher**
  reads' alleles per copy, not just the tied reads.

### What to build next (in order)

1. **Haplotagged teacher reads.** Re-call longcallR with BAM output, so each copy's teacher reads carry phase. This is light: about
   0.7 s per copy.
2. **The tier 1–3 screen as reported-only O3 outputs.**
   - Tier 1: bipartite-matching certificate.
   - Tier 2: novel recurrent haplotype.
   - Tier 2b: unexplained indel and tie widening by edit count.
   - Tier 3: collapse-compatible pair plus anchored ownership.
3. **Extend `donorsim.py` to simulate donor-only copies** (today it simulates only orthologs), together with absent copies and
   uneven allele expression.
4. **Pre-register the screen on a fresh HPRC donor.**
   - Metric: AUC, and sensitivity at the incumbent flag's own false-alarm rate (a matched operating point, no new constant).
   - Controls: EIF3C as a negative control (r1329); haploid X/Y copies in males as a free control (a het on a non-PAR Y copy is
     always a mixture).
   - HG002: reported, never claimed, because of its liftover-label artifacts.

### Related tonight

- **Naive O2 (`docs/O2_NAIVE_TIED_2026-10-22.md`):**
  - **Highest tied AS:** it is minimap2, adding nothing.
  - **Crediting every tied copy:** it recovers a few expressed copies that primaries miss, but invents more.
- **RNA-only O2 gates (`docs/O2O3_EXPLORE_2026-10-22.md`):** none brings the held-out donor under 2%. The hindsight H5 does, but at a
  forced cost (statement 2).


---

## O3_SCREEN_SPIKE_2026-10-22.md

*(original: `docs/archive/2026-10/merged/O3_SCREEN_SPIKE_2026-10-22.md`)*

EXPLORATORY (spike, post hoc)

## O3 RNA-only "look harder" screen, tiers 1-2: a 45-minute spike on HG002 and HG00733 (2026-10-22)

**What.** Tiers 1 and 2 of `docs/O3_FINGERPRINT_SURVEY_2026-10-22.md` (definitions: `docs/o3survey_2026-10-22/C_certificate.md` §1-3).
- Built as `bench/o3_graph/o3screen.py OUTDIR COPIES GENOME OUT.tsv`, one row per family (2 sanity tests: `test_o3screen.py`). It
  reads only the frozen O2-VC runs `/mnt/linuxdisk/tmp/o3_graph/o2vc/{hg002,hg00733}`: HG002 about 2.3 min and 220 MB, HG00733 4 s.
- The truth is used only for the family labels. Every variant is reported, including those made after looking (marked †).

### Definitions as built (constants: a fact seen by >= 2 reads, a conflict at >= 2 positions)
- **Informative positions of copy c:** every position of c in the family's `columns.tsv`, any kind, any tied read.
- **Pattern:** a read's bases at the informative positions it covers, base quality >= 10, in c's genome strand.
  - Teachers come from `teach/<family>.bam` (attributed to copies as `o2vc_data.scan` does); tied reads come from `columns.tsv`
    (the strand is told from the genome base against the row's `ref`).
  - A read covering no position has no pattern.
- **Conflict:** two patterns differ at >= 2 positions that both cover.
- **Classes: greedy first fit,** not exact match.
  - Distinct patterns are taken most positions first, then most reads; each joins the first class with no member it conflicts with.
  - A class of >= 2 reads is supported.
  - A class's pattern is its certified consensus: positions where >= 2 members carry the same base and it is the strict majority.
    The brief left this open.
- **Tier 1 (per copy):** a greedy multi-start clique of pairwise-conflicting supported teacher classes. It is a lower bound on
  the clique number. A copy is flagged when the clique is larger than its capacity.
  - Capacity: 1 on a non-PAR chrX/chrY copy of a male (PAR from `hg002truth.PAR`), else 2.
  - Sex was checked in the data. HG002 is male: tied reads in chrY's HSFY1; 25 of its 54 chrX/Y copies get capacity 1. HG00733
    is female: 0 of 127 non-PAR chrY copies have an ortholog, and 28 non-PAR chrX copies have both a mat and a pat one.
- **Tier 2 (per family):** the tied reads of `calls.tsv`, classed per anchor copy and counted once per read set.
  - A supported tied class is **novel** when it conflicts with every supported teacher class of every copy its reads tie to; each
    copy is compared in its own coordinates.
  - **t2 (literal):** a copy with no supported teacher class imposes nothing (vacuous).
  - **t2t (taught):** every candidate copy must have >= 1 supported teacher class.
- **Today's flag:** >= 1 row in `<family>.flags.tsv`.

**Variants** (flags of the script):
- **v1** is the default above.
- `--exclude-editing` drops `editing` rows, as C §1 drops A>G.
- **v2 `--pair-certified`:** the two conflicting positions must be carried *together* by >= 2 members of each class (C S1).
- **v3 `--merge-single`†:** classes one position apart are merged by single linkage (C §1). Added after diagnosing v1's cliques.
- **v4 = v3 + `--leak-closure`†:** tier 1 needs a Hall deficiency > 0 (C S4). Each clique class may also take the slots of any
  other copy with a teacher class it does not conflict with. Added after the HG00733 false alarms were traced to leaked teachers.
  - The comparison is made in the other copy's coordinates, through the `columns.tsv` homology.
  - A copy that no tied read links to c shares no position, so it counts as compatible (the conservative choice).
    `--leak-needs-overlap` (v4o†) requires >= 1 shared position instead.
  - A first leak-closure run compared positions without translating them, so it was vacuous. It was fixed before every number here.

### HG002 (94 families with tied reads; the 15 PAR families excluded; 73 "others")
Positives: (a) SPDYE1, RASA4B, USP17L11, POTEM, BMS1P1; (b) = (a) + NPIPB2.

| arm | (a) | (b) | others flagged (false alarms) |
|---|---|---|---|
| tier 1 v1 | 2/5 SPDYE1, RASA4B | 3/6 +NPIPB2 | 9/73: AC1275021, ANKRD36C, DUSP22, GTF2H2C, NBPF1, NOMO1, SMG1, SMN2, THOC3 |
| tier 1 v1, no editing | 1/5 SPDYE1 | 2/6 | 8/73 (SMG1 drops) |
| tier 1 v2 | 2/5 | 3/6 | 9/73 (same as v1) |
| tier 1 v3† | 1/5 SPDYE1 (3/2) | 2/6 +NPIPB2 | 5/73: AC1275021, ANKRD36C, GTF2H2C, NBPF1, SMG1 |
| tier 1 v1 + leak† (overlap required) | 1/5 SPDYE1 (1/5) | 1/6 (2/6 +NPIPB2) | 6/73: DUSP22, GTF2H2C, NOMO1, SMG1, SMN2, THOC3 (7/73 +NBPF1) |
| tier 1 v4† | 0/5 | 0/6 | 2/73: GTF2H2C (.1 deficiency 3), SMG1 (.2, 1) |
| tier 1 v4o† | 0/5 | 1/6 NPIPB2 (.9, deficiency 1) | 3/73: GTF2H2C, NBPF1, SMG1 |
| tier 2 v1 literal / taught | 0/5 / 0/5 | 1/6 NPIPB2 (31 reads) / 0/6 | 5/73 HSPA6, MSTRG3150, SMN2, THOC3, TRIM74 / 4/73 |
| tier 2 v2 literal / taught | 0/5 / 0/5 | 1/6 / 0/6 | 2/73 MSTRG3150, SMN2 / 1/73 SMN2 |
| tier 2 v3† (= v4) literal / taught | 1/5 RASA4B (119 reads) / 1/5 | 2/6 +NPIPB2 / 1/6 | 7/73 DUSP22, GOLGA8A, HSPA6, MSTRG3150, SMN2, THOC3, TRIM74 / 6/73 (no MSTRG3150) |
| today's flag | 1/5 SPDYE1 | 2/6 +NPIPB2 | 7/73: AC1275021, GOLGA8A, MSTRG3150, NCF1, NOMO1, SMG1, SMN2 |

The PAR families flagged (excluded from the table) were CD99 and SLC25A6 (tier 2) and SLC25A6 (today's flag).

**Ranking by tier-2 read support (v1, literal).**
- THOC3 2,887 > SMN2 1,594 > MSTRG3150 144 > HSPA6 109 > **NPIPB2 31** > TRIM74 20. Every other family scores 0.
- AUC-like score (positives against the 73 others, ties count 0.5):
  - tier-2 reads: 0.47 (a) / 0.55 (b);
  - today's flag reads: 0.56 (a) / 0.64 (b).
  - With nearly every family at 0, both are close to chance.
- **v3 ranking:** THOC3 2,892 > SMN2 1,602 > DUSP22 152 > GOLGA8A 146 > MSTRG3150 144 > **RASA4B 119** > HSPA6 109 > **NPIPB2 31** >
  TRIM74 20, then 0.
  - AUC-like score: tier-2 reads 0.55 (a) / 0.61 (b); today's flag 0.56 / 0.64.
  - SPDYE1 is caught only by tier 1 and today's flag.

**Why three positives are missed in every version.**
- USP17L11 and POTEM have **no taught copy**: 0 copies with a supported teacher class, from 12 and 26 teachers over 62 and 11 copies.
- BMS1P1 has 2 taught copies, 13 tied reads, and neither a clique > 2 nor a recurrent tied class.

### HG00733 (simulated, female, no donor-only copies: every flag is a false alarm; 174 families with tied reads)
The donor lacks a reference copy in 46 families, reproduced from `orthologs.tsv` against `pool_copies`.

| arm | false alarms / 174 | in the 46 absent-copy families |
|---|---|---|
| tier 1 v1, v1 no editing, v2, v3 | 3: GOLGA6L6, PRAMEF10, SPATA31A1 | 1/46 SPATA31A1 |
| tier 1 v1 + leak†, v4† | **0** | 0/46 |
| tier 1 v1 + leak†, v4o† (overlap required) | 1: GOLGA6L6 (.1, deficiency 1) | 0/46 |
| tier 2 v1 literal / taught | 10: FFAR3, KRT17P1, LIMS1, LYZL1, NCF1, OR1E1, PRAMEF11, RHOXF2B, TPSB2, ZC3H11A / 8 | 2/46 LIMS1, PRAMEF11 / 2/46 |
| tier 2 v1 no editing, literal / taught | 8 / 6 | 1/46 PRAMEF11 |
| tier 2 v2 literal / taught | 5: FFAR3, OR1E1, PRAMEF11, TPSB2, ZC3H11A / 5 | 1/46 PRAMEF11 |
| tier 2 v3† (= v4) literal / taught | 8 (v1's minus PRAMEF11, ZC3H11A) / 6 | 1/46 LIMS1 |
| today's flag | 5: FFAR3, GOLGA6L6, KRTAP21, TPSB2, ZC3H11A | 0/46 |

### Diagnosis (post hoc; it explains the variants)
- **v1's cliques are mostly fragmentation.**
  - THOC3.0 has 4,054 teachers and 508 informative positions. They give 109 supported classes: one of 3,685 reads and 108 of 2-9
    reads. The greedy clique is 33.
  - Of the 288 reads in the small classes, 267 differ from the main consensus at exactly 1 position.
  - So a large class sheds reads, because "non-conflicting with every member" adds two members' independent errors together.
    Each fragment then needs only one certified deviation to conflict with every other fragment.
  - v2 does not stop this, since each fragment carries both of its positions; v3 does.
- **HG00733's tier-1 false alarms are teacher leakage** (traced with `sim/truth.tsv`). The third class of every flagged copy is reads
  from a paralog that minimap2 placed as teachers:
  - GOLGA6L6.1: 2 reads from .3;
  - PRAMEF10.1 and .4: 5 and 2 reads from each other;
  - SPATA31A1.1 and .3: 2 and 4 reads from each other.
  Leak closure removes all three. GOLGA6L6 goes only through the vacuous rule: no tied read links .1 to .3. Requiring a shared
  position leaves it (v4o: 1/174).
- **Leak closure is per copy here.** Another copy's slots are not reserved for that copy's own haplotypes, so it is weaker than
  C's family-wide Hall test.
  - SPDYE1.18 is lost in v4 and in v4o: at least one of its 3 classes shares positions with a sibling copy's class without
    conflicting.
  - NPIPB2.9 survives only in v4o.
  - The family-wide matching, with the same haplotype seen from two copies counted once, was not built.
- **RASA4B.2** (v1) is two real haplotypes plus a 13-read class at two adjacent `editing` positions.

### Reading
1. As specified (v1), neither tier beats today's flag. Tier 1 catches 2/5 but false-alarms in 9/73 HG002 and 3/174 HG00733 families (bar: 0); tier 2 catches 0/5. Its cliques are fragmentation and teacher leakage, not haplotypes.
2. C's own rules, added post hoc† (merge, then leak closure), make tier 1 specific (v4: 0/174 HG00733, 2/73 HG002) but blind (0/5; v4o 1/6 via NPIPB2). Tier 2 v3 catches RASA4B (and NPIPB2 in (b)) but ranks CNV-prone SMN2, GOLGA8A and THOC3 above them (AUC-like 0.55 / 0.61 vs today's flag 0.56 / 0.64).
3. USP17L11, POTEM and BMS1P1 have no taught copy or too few reads, so no teacher-based tier can see them. Nothing here is claimable; next is a family-wide Hall test and a prereg on a fresh donor with simulated donor-only copies. Absent copies stay invisible (0/46).


---

## LANDING_LOCI_NEARBEST_2026-10-23.md

*(original: `docs/archive/2026-10/merged/LANDING_LOCI_NEARBEST_2026-10-23.md`)*

## Landing loci with near-best support, P2 control (2026-10-23)

**Question (user):** are there loci that only secondary alignments support, and are they copies the sample has ("phantom")?

**Data.** The port's landing-loci scan on the P2 control (53 real gorilla KB3781 families, unmasked reference; run
`/mnt/linuxdisk/tmp/rna_allele/p31/cand.loci.tsv`): loci outside a family's copies where >= 6 of its reads land (primary or
secondary). Annotation: `winloci_data/GGO_genomic.gff` genes / pseudogenes; other families' copies from `intervals.tsv`.

### Counts
| | loci |
|---|---|
| all landing loci | 746 in 40 families |
| only weaker (secondary, > 10 AS below the read's best) records | **722** |
| >= 6 reads within 10 AS of their best (near-best) | 13 |

**The 722 secondary-only loci** are reference loci the sample's reads resemble but do not come from: silent paralogs and pseudogenes,
as in the one-copy simulation (a silent pseudogene collected 516 secondaries, none near-best). They mark family members for O1; they
are not copies the sample expresses, and not O3 copies (a copy the reference lacks puts its reads on the nearest copy, primary or
tied, which the certified / putative tiers read).

### The 13 near-best loci
| class | loci | what they are |
|---|---|---|
| **another family's copy** | 7 (GWFAM169 x3, GWFAM401 x2, GWFAM402, GWFAM247) | the locus is a copy of GWFAM401 / 407 / 169 / 269: the reads are shared across two O1 families (an O1 boundary, the 16q NPIP-region families and TRIM5 / TRIM6 / TRIM34), not a new copy |
| **annotated gene, no family, strong** | 2: GWFAM123 at NC_073240.2:60.13-60.24 Mb (296 near-best of 516; MYO1E / FAM81A region incl. pseudogene LOC115930967); GWFAM99 at NC_073247.2:103,873,272-103,874,657 (151 of 323; PCDH11X with pseudogene LOC101144990) | expressed loci missing from the copy table: candidate unlisted copies (likely the annotated pseudogenes being transcribed) |
| **annotated gene, no family, weak** | 4: GWFAM169 / 401 at NC_073242.2:16.19-16.31 Mb (19 and 10), GWFAM268 (8 of 778), GWFAM28 (14 of 972) | a handful of near-best reads beside hundreds of secondaries: mostly paralog spill-over, a few reads possibly from the locus |
| **unannotated** | 0 | |

**Answer.** Secondary-only loci exist in large numbers (722 / 746) and are silent reference paralogs / pseudogenes, not phantom copies.
Of the 13 with near-best support, none is unannotated: 7 are other families' copies (O1 boundary), 2 are strong candidates for
expressed copies missing from the copy table (both overlap annotated pseudogenes), 4 are weak. Nothing here is a copy absent from
the reference. Not followed up: whether the GWFAM123 / GWFAM99 loci belong in those families' copy tables (O1 completeness).


---

## O3_FIRST_PILOT_2026-10-23.md

*(original: `docs/archive/2026-10/merged/O3_FIRST_PILOT_2026-10-23.md`)*

EXPLORATORY (pilot, post hoc; no pre-registration, nothing claimed)

## O1 -> O3 -> O2 pilot on HG002, with Strainy and the advisor's silent-copy check (2026-10-22/23)

**Code:** `bench/o3_graph/o3first_pilot.py` (`copies`, `silent`, `strainy-inputs`, `score`). **Work dir:** `/mnt/linuxdisk/tmp/o3_graph/o3first`.
**Inputs:** HG002 Kinnex reads on CHM13 (`hg002proof/chm13.bam`); the O2-VC copy table (`o2psv_dev/pool_copies.tsv`). The table lists
every exon variant of a copy (793 of its 1,686 copies overlap), so each copy is merged into its exon union; chrM copies are dropped.
Families: the 94 with tied reads from the 2026-10-22 spike, with the same labels: positives (a) SPDYE1, RASA4B, USP17L11, POTEM,
BMS1P1; (b) = (a) + NPIPB2; 73 others; the 15 PAR families are excluded.

### 1. O3 before O2: family-level detection (every arm on the same 79 families)
- `candidate_copies` is the shipped opt-in O3 stage, run with its defaults on each family's whole read net before any tied read
  is placed: 7 batches, 4-6.5 min each, 14.4 GB peak.
- Strainy 1.2 (bioconda) ran as `-s phase -m hifi`, with the family's copies (spliced exon unions) as the collapsed
  sequences and every read with a record on them. `--max-unitig-coverage` was raised to 100000 because RNA depth exceeds the
  500 default; every other setting was left at Strainy's default, including 20x minimum coverage, 1 kb minimum length and MAPQ >= 20.
  - A copy's haplotypes = the most Strainy clusters (>= 2 reads each, unclustered excluded) that cover one of Strainy's SNP
    positions together.
  - A copy is flagged when that count is larger than its allele slots: 2, or 1 for a non-PAR X/Y copy, since HG002 is male.

| arm | (a) /5 | (b) /6 | others flagged /73 |
|---|---|---|---|
| Strainy: haplotypes > allele slots | 1 SPDYE1 | 2 (+NPIPB2) | 7: NBPF1, AC1275021, CBWD2, SMG1, GTF2H2C, ANKRD36C, DUSP22 |
| candidate_copies: flagged candidate | 1 SPDYE1 | 1 | 14 |
| candidate_copies: any candidate | 4 SPDYE1, RASA4B, POTEM, BMS1P1 | 4 | 20 |
| spike tier 1 (10-22) | 2 SPDYE1, RASA4B | 3 | 9 |
| today's flag (O2-VC, 10-22) | 1 SPDYE1 | 2 | 7 |

- **Strainy coverage.** Strainy phased >= 1 copy in 50 of the 79 families. USP17L11 and POTEM have no copy at >= 20 reads, so
  Strainy never phased them. GOLGA8A hit the 280 s limit. SPDYE1.20 and NPIPB2.16 each carry 3 concurrent haplotypes.
  RASA4B's 2 phased copies have 1 each, and BMS1P1's 3 phased copies have at most 2.
- **Labels.** 3 of Strainy's 7 "false alarms" (CBWD2, NBPF1, GTF2H2C) are families where the truth labels hundreds of reads as
  coming from donor copies with no CHM13 counterpart: 395, 243 and 212 reads. They were not among the verified positives and may
  be liftover artifacts. Every arm's count moves if they are real, so the label set must be verified before any arm is compared.
- **Candidate support.** candidate_copies' unflagged candidates in RASA4B, POTEM and BMS1P1 have 3-4 reads each, the same scale as
  its candidates in the other families.

### 2. The advisor's claim: a silent copy is found through its family's multimapping reads
- **Definitions.** A silent copy has no read with its primary record on it; an expressed copy has >= 3. A silent copy is
  reached when a read of an expressed copy of its family has a secondary record with >= 100 bp on it. Identity and coverage
  come from minimap2 -c -x asm20 against the nearest expressed copy.
- **Result.** 236 silent copies in 75 families (NUMT families excluded); 179 (76%) are reached, 166 (70%) by >= 2 reads.

| identity | silent | reached |
|---|---|---|
| >= .99 | 66 | 66 (7 of them fragments covering < 50% of the gene) |
| .97-.99 | 41 | 32 |
| .90-.97 | 104 | 64 |
| < .90, or no asm20 hit | 25 | 17 |

Expression is not the limit, divergence is. These are annotated copies; unannotated ones were not measured.

### Not done
- The read-level O2 effect: whether the donor-only copies' reads move to the candidates, and the false moves.
- Verifying the labels of CBWD2, NBPF1, GTF2H2C and WASHC1.
- A Strainy run with a lower coverage floor.

### 3. Label check (2026-10-23): which `no_counterpart` reads come from sequence CHM13 lacks
- **Method** (`o3first_pilot.py verify-queries` / `verify-call`, then the per-base pass).
  1. Each `no_counterpart` read takes its best donor alignment (mat/pat, tied and untied).
  2. Reads within 10 kb of each other form donor loci; each locus with >= 2 reads becomes a query, the locus plus 20 kb of flank
     on each side from hg002v1.1.
  3. The queries are aligned to CHM13 with minimap2 asm20 `--eqx`.
  4. A read comes from donor-only sequence when >= 50 of its exon bases fall in an insertion (>= 50 bp) of the alignment that
     spans both flanks.
  - A first version that asked only whether the locus aligned across both flanks called every locus "counterpart". asm20
    bridges multi-kb insertions inside one gapped alignment, so that version was discarded.

| family | finding | labelled reads from donor-only sequence |
|---|---|---|
| SPDYE1 | 6.4 kb donor insertion (chr7_MATERNAL 104,183,399-104,189,258) | 99 / 342 |
| POTEM | 0.35 kb insertion | 2 / 12 |
| RASA4B, SDHA | 6.3 kb and 0.9 kb insertions | 1 / 19, 1 / 10 |
| NPIPB2 | the known 57-bp insertion in a copy CHM13 has (allelic structural variant) | 25 / 43 |
| NBPF1 | a 6.4 kb donor DELETION: HG002 lacks a unit CHM13 has; 11 reads span the junction | 0 (11 junction reads) |
| CBWD2, GTF2H2C, WASHC1, USP17L11, BMS1P1, AGAP4, DUXAP9, FRG1, CTAGE15, ... | allele-level counterpart (<= 1% divergence, no insertion): liftover artifacts | 0 |

**HG002 has one solid extra-copy positive, SPDYE1. Its 2026-10-22 positive set was partly wrong:** USP17L11 and BMS1P1 have no
read from donor-only sequence, and RASA4B has one. The family-level tables above cannot rank arms.

### 4. Read level: does any arm find SPDYE1's hidden copy?
- **minimap2.** All 99 reads have their primary on CHM13 copy SPDYE1.18, unique, among ~1,050 reads there (9%). Their median NM
  is 7 against CHM13 and 2 against the donor.
- **The hidden copy's sites.** At the columns where its reads all carry a non-reference base, the donor's own SPDYE1.18 reads
  mostly carry it too (alleles shared with one of its haplotypes). It is separable only by linkage across columns.
- **candidate_copies.** 37 of the 99 reads sit in clusters linked to a reference copy as alleles (delta 0.96%), and 62 in
  in-reference clusters. Its flagged SPDYE1 candidates (81 + 11 reads) contain none of them.
- **Strainy.** 98 of the 99 fall in SPDYE1.18's main group (1,018 reads). Its pileup caller's allele-frequency floor (0.2) is
  above the hidden copy's 9%. Its SPDYE1 flag came from copy SPDYE1.20.

**Every arm "caught" SPDYE1 for another reason.** What the case needs:
- variant columns admitted by recurrence or significance, not by an allele-frequency floor (IsoCon's statistical test);
- reads linked across columns with zero disagreement and refined per cluster (Strainy's graph);
- the cluster's consensus as the new copy (IsoCon).


---

## O3_ONE_COPY_SIM_2026-10-23.md

*(original: `docs/archive/2026-10/merged/O3_ONE_COPY_SIM_2026-10-23.md`)*

## O3 on one erased copy, simulated (2026-10-23)

**Why:**
- The user asked for one case instead of a 60-family panel: "can we somehow only focus on one erased copy? We could even simulate it
  so it's clean".
- They then named the base case: "absent copy that is hidden in multi mapping reads ... [the advisor] believes multi mapping reads
  are everywhere there is similarity so pseudogenes should be trackable".
- They chose to run both cases: the copy missing from the genome first, then the copy present but unannotated.

**Code:** `bench/o3_graph/o3sim_one.py` (pick / sim / report). Work dir: `/mnt/linuxdisk/tmp/o3sim1`.

### Setup
- **Family:** gorilla KB3781 GWFAM248, 5 copies, with GWFAM248:3 erased. The selection rule was fixed before any port run:
  - the simplest family (fewest copies, >= 3) whose erased copy's transcript is >= 97% identical over >= 80% of its length to
    >= 2 survivors;
  - >= half of the erased copy's reads must be multimappers, i.e. have a record within 10 AS of their best on >= 2 survivors.
  - The erased copy is 97.4% identical to survivors :1 and :2. 85% of its pick reads tie on >= 2 survivors.
- **Truth transcripts:** each copy's dominant fibroblast isoform (the most frequent intron chain of the fibroblast BAM's primary
  reads), taken on both haplotypes (A = `_pri`; B = lifted to the other haplotype, `isocon_sim.lift_blocks`).
  - In KB3781, :3 is homozygous (A = B).
  - :1 and :2 differ at 1 and 0 positions.
- **Reads:** 100 per allele per copy. The model is `isocon_sim`'s: 5' / 3' truncation, 0.5% substitutions, 0.25% insertions,
  0.25% deletions, seed 1. Every read's origin is known.
- **Reference:** a mini reference, i.e. the family's loci +- 50 kb, aligned with the Amendment 7 command (`-N 50 -p 0.1`).
- **Port:** v5 defaults, on the copy table without :3.

### Results
| case | erased copy's reads | the port | where the family's reads land outside the copy table |
|---|---|---|---|
| **missing** (masked in the genome, expressed) | 200 of 200 on survivors; 163 are multimappers (>= 2 piles) | one haplotype of 184 reads, **all from :3**; it takes the free second-allele slot of survivor :1 and is **not certified** | an unlisted similar locus, 473 secondary-only reads |
| missing, survivors heterozygous (**synthetic**: 6 substitutions on :1 / :2's B alleles) | the same | the same 184-read haplotype is **certified** (35 private alleles; 0 survivor reads in it); contig 1,306 bp = :3's transcript (1,623 bp) at **99.92%** (0 mismatches, 1 indel) | the same |
| **unannotated**, expressed (in the genome, not in the copy table) | its own locus (200 reads within 10 AS of their best) | nothing to certify on the survivors | **on :3: 715 reads** (200 near-best, 515 secondary) |
| unannotated, **silent** (a pseudogene: no reads of its own) | none | nothing | **on :3: 516 secondary-only reads** from the survivors' reads |

### What it shows
1. **The print of the missing copy is there and the port finds it.**
   - The missing copy's reads hide as multimappers across the survivors, as the advisor predicts.
   - Clustering by allele pattern pulls all of them into one haplotype that contains nothing else.
   - On survivor :1 that haplotype carries a non-reference allele at 26 of 27 informative columns; :1's own haplotype carries 0.
2. **Proving it is extra depends on the survivors (Theorem A).**
   - Where the survivors' alleles are distinguishable, both of their slots are full and the haplotype is certified, with a 99.9% contig.
   - KB3781's real survivors are homozygous here. A homozygous survivor has a free slot, and RNA cannot rule out that the haplotype
     is that survivor's second allele, so it is not certified.
   - That case is a putative copy: a slot-filler that disagrees with the copy's reference almost everywhere.
3. **Pseudogenes in the genome are trackable from multimappers alone.** The advisor's claim holds here:
   - Silent :3 collects 516 secondary alignments from its paralogs' reads.
   - So does an unlisted similar locus near the family (280-473 reads).
   - The signal is a pile of secondary alignments at a locus outside the copy table; none of them is within 10 AS of the read's best.
4. **A copy missing from the genome and silent leaves no reads,** so RNA has nothing to find (not simulated).

### Reproduce
```
O=/mnt/linuxdisk/tmp/o3sim1
python3 bench/o3_graph/o3sim_one.py pick --out $O
python3 bench/o3_graph/o3sim_one.py sim --out $O --name missing --case missing
python3 bench/o3_graph/o3sim_one.py sim --out $O --name missing_shet --case missing --s-het 6
python3 bench/o3_graph/o3sim_one.py sim --out $O --name unannot_expressed --case unannotated
python3 bench/o3_graph/o3sim_one.py sim --out $O --name unannot_silent --case unannotated --e-reads-per-allele 0
python3 bench/o3_graph/o3sim_one.py report --out $O --name <name>
```

### Sweep (port v6; `o3sim_one.py sweep`, 40 runs, `/mnt/linuxdisk/tmp/o3sim1/sweep.tsv`)
**Setup.**
- **Another individual** (the user: "the method should somehow account for individual differences of not having used the same
  individual for the genome and for the RNA samples"): the sample's own variants, which the reference lacks. Per copy, 0.1% of
  sites change on both alleles and 0.1% on one allele, a human-like rate. This is the CHM13-vs-HG002 situation.
- **Divergence:** E's 34 differences from its closest survivor :1 are kept with probability f, the rest reverted to :1's base.
  f = 0.5 / 0.25 / 0.1 / 0.05 gives about 1.3 / 0.65 / 0.26 / 0.13% divergence.
- **Heterozygous survivors** (synthetic): survivors whose real alleles differ at < 2 positions get a B allele with 6 substitutions.

| E (the erased copy) | same individual, survivors as in KB3781 | same individual, survivors heterozygous | other individual, survivors as in KB3781 | other individual, survivors heterozygous |
|---|---|---|---|---|
| 100 / 30 / 10 / 5 reads per allele | **putative** (pure print; contig 99.9-100%) | **certified** | certified at 100 / 30 / 10; missed at 5 | certified at 100 / 30 / 10; missed at 5 |
| 3 reads per allele | missed (below the 6-read floor) | missed | missed | missed |
| f 0.5 / 0.25 (1.3 / 0.65%) | silent in :1's free slot | certified | certified | certified |
| f 0.1 (0.26%) | silent | count right, but :1's own allele is named extra | certified | count right, wrong one named |
| f 0.05 (0.13%) | silent | count right, wrong one named | silent | count right, wrong one named |
| control (nothing erased) | 0 false calls | 0 false calls | 0 false calls | 0 false calls |

**What the sweep shows.**
1. **The reference being another individual does not create false calls (0 in every control).** The certificate never compares
   with the reference sequence to decide:
   - a column is informative only when the sample's own reads carry >= 2 alleles, so a variant the sample carries on both
     alleles is invisible to it;
   - haplotypes, the count against slots, and the putative-copy flag (tie patterns) all compare the sample with itself.
   Only two parts read the reference:
   - the slot order (closest to the reference first) decides which haplotype is named extra when the count exceeds the slots;
   - the DADA2 rule's reference parents.
2. **Another individual helps more than it hurts here.** The sample's own heterozygous sites make KB3781's homozygous survivors
   heterozygous, which fills their slots, so the erased copy becomes certifiable. It costs depth, though: the erased copy's two
   alleles now differ, so each needs 6 reads on its own (missed at 5 per allele).
3. **The limits are depth and divergence.**
   - Depth: about 5-10 reads per allele.
   - Divergence: below about 0.3% from its closest survivor (<= 4 differences), the extra copy is still counted, but RNA cannot
     say which of the three haplotypes is the extra one, so a survivor's allele can be named instead.
4. **The sweep found a real bug (v6).** Before the fix, the conflicting set counted 2-3-read artifact clusters. Three of them beat
   a real 93-read allele, freeing its slot for the erased copy (`expr100_6_other`: putative, not certified). The 6-read floor
   now applies before the count. On real development data this fix moved P1 from 5,762 to 6,129 D right (passes) with P2 still
   0 / 53 (`docs/PREREG_o3_port_2026-10-23.md`, development log).

### Correlation clustering in place of label propagation (`--clusterer cc`; sweep in `/mnt/linuxdisk/tmp/o3sim1_cc`)
Same 40 runs; the cells that change:

| run | label propagation | correlation clustering |
|---|---|---|
| other individual, 100 reads per allele: the erased copy's pure haplotype | 110 reads | **179** |
| other individual, 30 / 10 per allele | 22 / 12 reads | 49 / 16 |
| other individual, 5 per allele | missed | **certified** (7 reads) |
| same individual, 10 / 5 per allele | 15 / 6 reads | 19 / 9 |
| controls (nothing erased) | 0 false calls | 0 false calls |
| erased copy <= 0.26% from a survivor | count right, wrong one named | the same |

**Why it helps.** When the sample is another individual, the erased copy's two alleles differ at a site or two.
- Label propagation links only reads that disagree nowhere, so it splits the copy into two half-size haplotypes.
- Correlation clustering keeps them together (one differing column is one observation).

The divergence limit is the slot logic's, not the clustering's.


---

## O3_PORT_RESULTS_2026-10-23.md

*(original: `docs/archive/2026-10/merged/O3_PORT_RESULTS_2026-10-23.md`)*

## O3 port (Strainy's phasing + IsoCon's statistics): results on gorilla KB3781 (2026-10-23)

**Status: pre-registered** (`docs/PREREG_o3_port_2026-10-23.md`, ddc1b4b7). The code was frozen for P3 at 1b4535cc.
**Verdict:** P1 fails, P2 passes and P3 H1 passes, so the port is **not adopted** as the O3 stage: "detects but not sensitive"
(H3 holds). Design: `docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md`. Code: `bench/o3_graph/o3port.py`
(11 tests).

### What the frozen port does
- **Strainy side**, per reference copy:
  1. **Piles:** a read joins every copy it aligns to within 10 AS of its best.
  2. **Informative columns:** >= 2 alleles with >= 2 carriers each, the minor allele significant under the reads' own empirical
     error rates (IsoCon), Bonferroni, alpha 0.01.
  3. **Clusters:** a zero-disagreement read graph, label propagation and per-cluster refinement.
  4. **Haplotypes:** the largest set of haplotypes that differ pairwise at >= 2 pile-level columns.
  5. **Chimeras and duplicates** (DADA2's rule, before slots): a haplotype that is a left/right join of two parents takes no
     slot, and one equal to a more abundant haplotype shares that haplotype's slot. Parents are every family copy's reference
     and the more abundant haplotypes, mapped across copies by asm20 homology. A haplotype equal to a copy's reference is that
     copy's allele.
  6. **Certificate:** allele slots (2 per copy) filled reference-like first by b-matching. A haplotype is certified when it gets
     no slot, carries >= 2 private alleles and has >= 6 reads (the shipped stage's own flag floor).
- **IsoCon side:** the shipped stage's reference-free candidates, kept when most of their reads reached the net without
  aligning to any copy.

### Results (the same reads and the same scorer as the shipped stage K)
| rule | port | K (`candidate_copies`, Amendment 15) | bar | port verdict |
|---|---|---|---|---|
| P1 deletion panel (53 families, dev): D right | 3,995 | 7,882 | >= 5,842 | **FAILS** |
| P1: false moves | 7 / 41,727 = 0.02% | 0.77% | <= 5% | holds |
| P2 no-deletion control (53, dev): families with a false flag | 0 | 10 | <= 8 | **PASSES** |
| P3 held out (30 families, frozen code): H1 families with a false flag | 3 | 6 | <= 4 | **PASSES** |
| P3: H3 false moves | 178 / 27,757 = 0.64% | 0.26% | <= 5% | holds |

- **Held-out flags.** The port flagged 6 haplotypes in 4 families:
  - 1 is a real reference-absent copy: GWFAM390, 8 reads, 100% identical to a maternal-only locus (CM054598.2:30.25 Mb).
  - 3 are B-haplotype alleles of existing copies: GWFAM205 and GWFAM256, at 99.94-99.96%.
  - 2 are unmatched: GWFAM256 (6 reads) and GWFAM318 (IsoCon side).
- **On P1.** The port gives back half the stage's deleted-copy reads, with about 38x fewer false moves.
- **On H3.** 362 held-out reads land on their own candidate. Many of the 178 false moves are reads of the existing copies' B
  alleles, which the b_allele contigs attract.

### How v1 became v4 (development on P1 / P2 only; every change uses an existing constant)
| version | change | P2 families flagged | P1 families with a D candidate |
|---|---|---|---|
| v1 | base design | 35 (97 flags) | 16 |
| v1 + empirical error | IsoCon's per-read error rates, not PacBio QVs (capped at Q93) | | |
| v1.1 | conflicts proven only on pile-level columns | 23 | |
| v1.2 | >= 2 private alleles | 22 (C1' classified: 45 c, 5 b, 2 pri, 1 a) | |
| v2 | family-wide chimera rule before the slots | 19 | 13 |
| v3 | the stage's flag floor (6 reads); slot order by divergence (support order: 396 S reads certified on P1) | 0 | 8 |
| v4 | + IsoCon side (via-B reads) | 0 | 16 (P1 harness: D right 3,995) |

### What we learned
1. **Full-length reads were not the problem.** 68 of 72 surplus haplotypes were made of near-full-length reads. They were not
   A-to-I editing either: their alleles are transition-rich and none was all A>G.
2. **RNA piles are full of small artifact clusters** (2-7 reads): chimeric molecules, some joining paralogs, and reads that fit
   poorly. They match neither parent (median identity 0.991). Without a support floor, every allele-slot certificate drowns in
   them.
3. **The certificate is honest but conservative (Theorem A in practice).** A haplotype whose reads also tie to a copy with a free
   slot cannot be proven extra. That is how GWFAM175's real paternal-only copy (181 reads) was lost on development. On held-out
   data the port still found GWFAM390's maternal-only copy, with 1 extra false family fewer than the bar allows.
4. **Read support and divergence answer different questions.** Filling slots by support certifies the survivors' own minor
   alleles; filling them reference-like first certifies the collapsed copy.
5. **Sensitivity is lost in two places.**
   - Deleted copies whose reads tie to a survivor with a free slot.
   - The 31 of 53 deleted copies whose reads never touch a survivor. Those depend on the IsoCon side, and that side is only as
     sensitive as the stage's reference-free clusters.

### Reported, not decided (the pre-registration's list; computed after the verdict, verdict unchanged)
Source: `o3port.py run --dump` (internal states as JSON lines; the products do not change: P2 and the P3 flagged families
reproduce byte for byte). Work dir `/mnt/linuxdisk/tmp/rna_allele/p16`.
- **Runtime** (one Python process): P1 193 s / 1.74 GB peak, P2 222 s / 1.80 GB (53 families each); HG002 SPDYE1 14 s on the
  stage's 1,000-read net, 47 s / 1.07 GB on its whole 3,200-read net.
- **Editing.** No certified flag depends on an A>G / T>C column: without them every certified haplotype still has >= 2 private
  alleles (P1 14 of 14, P3 5 of 5 Strainy-side flags). Editing-type alleles are 9-40% of the private alleles, as expected for transitions.
- **Per-copy haplotype count in the control (P2) against mat / pat** (`control/copies_lift.tsv`; slot holders the port assigns
  per copy, all 201 copies have a pile):

  | truth | 0 holders | 1 | 2 |
  |---|---|---|---|
  | T2d (alleles differ), 174 copies | 27 | 47 | 100 |
  | T2i (alleles identical), 24 | 13 | 10 | 1 |
  | T1-candidate (one haplotype only), 1 | 0 | 1 | 0 |
  | T? , 2 | 0 | 1 | 1 |

  The port separates the two alleles of 100 of 174 heterozygous copies (median 35 exonic differences between them; merged ones
  median 19) and over-splits 1 of 24 homozygous copies.
- **HG002 SPDYE1 (development case, human).** The hidden copy's 99 reads (pilot §4) are all in the stage's net, but only 37 are in
  its 1,000-read sample. **The port does not certify the hidden copy:**
  - on the stage's net, 13 of the 37 sit in a 29-read cluster outside the conflicting set and 5 in a 22-read haplotype the DADA2
    rule calls a chimera;
  - on the whole net (report-only), 39 of the 99 form a 94-read haplotype the rule calls a join of another haplotype and
    SPDYE1.18's reference: the copy is a mosaic of the donor's own alleles (pilot §4: "separable only by linkage across
    columns"), which is exactly what a chimera looks like. A different 99-read haplotype on SPDYE1.20 is certified (0
    insertion reads; unverified).
- **Reproducibility fix.** The consensus broke 50/50 ties by Python's per-process string-hash order (2 of 29 P1 contigs differed
  by 1-2 bases between runs; read sets and certificates unaffected). Reads are now taken in sorted order; checked across hash
  seeds.

### Where the port loses deleted-copy reads (P1 development data, post hoc)
The furthest stage each of the 6,414 used D reads reaches, and how many of them sit in the shipped stage K's flagged candidates:

| stage | D reads | in K's flags |
|---|---|---|
| no pile (best alignment > 10 AS off every copy, or no record on a copy) | 2,265 | 2,154 |
| in a pile, no cluster (cover no informative column) | 1,222 | 433 |
| cluster outside the conflicting set | 399 | 361 |
| merged away / below 6 reads / duplicate / chimera | 154 | 115 |
| took a free allele slot (Theorem A) | 1,343 | 1,190 |
| certified | 1,031 | 988 |

- **No pile.** 976 of the 2,265 have no record on a copy (via B); the IsoCon side holds 752 of them. The other 1,289 align
  better elsewhere. The IsoCon side's v4 rule (members >= half via B) does not reach them; a rule on the actual piles (members
  >= half in no pile) would hold 2,109 D reads on P1 (v4: 804) with the same 5 S reads, and adds no flag on P2.
- **Free slots.** 10 deleted copies' haplotypes took a slot on a survivor, 9 of them on survivors whose own alleles differ
  (T2d) but were not separated. Divergence does not rescue them: the deleted-copy haplotype differs from the survivor's at
  a median of 27 informative columns (3-232), but true allele pairs reach 36 (P2) and 99 (P1 survivors), and some copies'
  mat / pat alleles differ at 150-240 exonic sites. RNA cannot tell an allele from an extra copy here (Theorem A).
- **No cluster.** In GWFAM401 all 275 D reads cover only columns no survivor read covers: no variation to phase against.


---

## PREREG_o3_port_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_port_2026-10-23.md`)*

## Pre-registration: O3 port (Strainy's phasing + IsoCon's test), gorilla KB3781 (2026-10-23)

Committed before the port has been run on any panel. Design: `docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md`.
The bars are the shipped stage's own bars (`candidate_copies`, Amendments 13-15: register rows 1221-1231 and 1247), so that the
port can replace the stage on equal terms.

### Arms, the same reads and the same scorer
- **P:** the port, run as designed.
- **K:** `candidate_copies` as it ran under Amendment 15. Its products are already on disk: a13a15 / a14a15 / a15h. It is not
  re-run.

### Decision rules
- **P1, development (deletion panel, a13a15).**
  - **Pass:** D right >= 0.80 x C = 5,842 (C = 7,303, Amendment 13b's comparator) AND false moves <= 5% of S reads.
  - **K on the same panel:** 7,882 and 0.77%, which passes.
- **P2, development (no-deletion control, a14a15).**
  - **Pass:** families with >= 1 false flag (class b / c / pri against mat / pat, Amendment 9's rule) <= 8 of 53.
  - **K:** 10 of 53, which fails.
- **P3, held out (a15h; run ONCE, after P1 and P2, with frozen code).**
  - **Pass:** families with >= 1 false flag <= 4 of 30 (H1) AND false moves <= 5% (H3).
  - **K:** 6 of 30 and 0.26%; H1 fails.
- **Adopt the port as the O3 stage iff P1, P2 and P3 all pass.**
  - P passes P1 but fails P2: "detects but not specific", not adopted.
  - P fails P1: "not sensitive", not adopted.

### Rules for development
- P1 and P2 are development. A change after looking is allowed only if it adds no new constant, is recorded in the doc with its
  reason, and is re-run on BOTH P1 and P2.
- Then the code is frozen (commit hash written here) and P3 is run once.

### Reported, not decided
- The per-copy haplotype count in the control against mat / pat: allele-level truth for every copy.
- How many certified flags depend on an editing-type column.
- Runtime.
- HG002 SPDYE1 (development case, human).

### Development log (P1 / P2 are development; every change below adds no new constant)
- **v1** (1cb54fa7), on GWFAM175 / GWFAM4 only:
  - per-read empirical error rates (NM / aligned bases, IsoCon) instead of PacBio QVs, which are capped at Q93 and let 2 carriers
    of a recurrent error pass;
  - a conflict needs >= 2 differing columns.
- **v1 deletion panel** (53 families), descriptive: 91 certified haplotypes in 35 families. Of their clusters, 40 hold only D
  reads, 47 only S reads and 4 both. 16 of the 53 deleted copies are caught: in 31 of the 53, no D read aligns to a survivor
  (A11).
- **v1.1:** conflicts are proven only on pile-level columns. Columns that pass only inside a refined cluster can come from
  recurrent errors or chimeras.
- **v1.2:** an extra copy needs >= 2 private alleles relative to the slot holders (mosaics excluded).
- **P2 on v1.2** (full classification, `ACC=a14 A=p14`): 53 flags in 22 of 53 families.
  - Classes: c_unmatched 45 (mostly 2-7 reads), b_allele 5 (219 / 118 / 97 / 37 / 3 reads), pri 2, a_haplotype_only 1
    (GWFAM175, 181 reads).
  - **C1' 22/53, FAILS** (bar <= 8). The shipped stage scored 10/53.
- **Diagnosis.** Of the surplus haplotypes, 4 of 72 are made of short reads (median read length 0.79-0.85 of the pile's p90),
  so requiring full-length reads would not remove them. Their private alleles are transition-rich and none are all A>G, so they
  are not editing. 28 of the 72 are mosaics of the copy's two alleles.
- **v2:** the family-wide DADA2 rule runs before the slots. A haplotype that equals a more abundant haplotype, or is a left/right
  join of two parents, takes no slot. Its possible parents are every copy's reference and the more abundant haplotypes, mapped
  by asm20 homology. Result: P2 19/53, P1 13 families with D.
- **v3:** MIN_SUPPORT 6, the shipped stage's own flag floor (candidate_copies --min-support 6, ruling R1), so the comparison is
  on equal terms: a haplotype needs 6 reads to take a slot or be certified.
  - Slot order compared on development: by divergence, P1 has 14 clusters, 1,031 D reads and 2 S reads; by read support, 301 D
    and 396 S. Divergence was kept.
  - P2 had 0 flags under both orders. GWFAM175's real paternal-only copy is no longer certified: its reads tie to a copy with a
    free slot, which is Theorem A's ambiguity.
- **v4 = FROZEN (1b4535cc):** v3 plus the IsoCon side, i.e. the stage's flagged candidates whose members are >= half via-B (unaligned)
  reads. On P1 that adds 13 candidates in 12 families (804 D, 5 S reads); on P2, 0.
- **P1 (v4, harness ACC=a13 A=p13): D right 3,995 (54.7% of C) < 5,842: FAILS.** False moves 7 / 41,727 = 0.02% (K: 7,882 and
  0.77%).
- **P2 (v4): 0 flags in 0 / 53 families: PASSES** (K: 10 / 53).
- **Verdict so far: "detects but not sensitive", not adopted** (pre-registered reading of a P1 failure). P3 is run once on the
  frozen code to record held-out specificity (H1 / H3); it cannot change the verdict.
- **P3 (frozen 1b4535cc, run once, harness ACC=a15h A=p15h):** 6 flags in 4 / 30 families.
  - **H1 3 / 30, PASSES** (<= 4). Classes: 1 a_haplotype_only (GWFAM390, 8 reads, maternal-only locus, identity 1.0000),
    3 b_allele, 2 c_unmatched.
  - **H3 178 / 27,757 = 0.64%, HOLDS** (<= 5%). 362 reads go to their own candidate.
- **FINAL: P1 fails, P2 and P3 pass. Not adopted ("detects but not sensitive").** Register row added.

### Development after the verdict (v5+, gorilla fibroblasts only; P1 / P2 are development, P3 is spent)
- **Loss accounting** (P1, `--dump`, results doc): of 6,414 used D reads, 2,265 sit in no pile (1,289 of them have a record on a
  copy but align better elsewhere, out of the v4 IsoCon side's reach), 1,343 took a free allele slot (Theorem A), 1,222 cover no
  informative column.
- **v5:** the IsoCon side keeps the stage's flagged candidates whose members are >= half in no phased pile (the reads the
  Strainy side never sees), not >= half via B. No new constant. `--isocon-rule via-b` reproduces v4 exactly (checked on P1).
  - P1 (harness ACC=a13 A=p17): D right **5,762** (78.9% of C) < 5,842: **FAILS by 80 reads**; false moves 7 / 41,727 = 0.02%.
  - P2: 0 flags in 0 / 53: passes.
- **The held-out for any version after v4 is a NEW gorilla fibroblast panel** (user 2026-10-23: "we need to stick to fibroblasts
  in gorilla"; no human test bed), registered before it is built.
- **v5.1** (`--absorb`, opt-in, **rejected**): design step 5 as written: a cluster outside the conflicting set joins the one kept
  haplotype it shares a column with and does not conflict with, and a haplotype merged across piles keeps every pile's reads.
  - P1: flags hold +87 D and **+140 S** reads (a 118-read survivor haplotype certified on GWFAM269); harness D right **5,825**
    (79.8% of C) < 5,842: fails by 17; false moves 20 (0.05%).
  - P2: **4 flags in 4 / 53** families (v5: 0).
  - Precision for 63 reads; v5 stays the candidate. Development stops here: neither version clears P1, and tuning further to pass
    it would be fitting the bar.
- **Putative copies** (2026-10-23, report-only beside the certificate; first called "look harder", renamed at the user's request): a slot holder >= half of whose reads tie to another copy as well
  as to the slot's, while < half of the reads of the slot copy's other holders do. It compares the sample with itself (tie
  patterns), never with the reference sequence. Products `PREFIX.putative.tsv` / `.fa`; and `PREFIX.loci.tsv` (register 1377).
- **v6** (found by the one-copy simulation with another individual's variants, `docs/O3_ONE_COPY_SIM_2026-10-23.md`): only clusters
  at the support floor (6 reads, MIN_SUPPORT) enter the conflicting set. Before, the largest set of pairwise-conflicting clusters
  counted 2-3-read artifacts too, and three of them beat a real 93-read allele, freeing its slot for the erased copy. No new
  constant: the floor's own meaning ("a haplotype needs 6 reads to take a slot"), applied before the count. `--v5-clique` restores v5.
  - P1 (harness ACC=a13 A=p21): D right **6,129** (83.9% of C) >= 5,842: **PASSES**; false moves 7 / 41,727 = 0.02%.
  - P2: **0 flags in 0 / 53**: passes. Putative: 2 families, one of them GWFAM175's real paternal-only copy (181 reads; lost to
    a free slot since v3), the other a 6-read flag (GWFAM4).
  - v6 passes both development panels. It needs a held-out (fresh gorilla fibroblast families) before any claim.
- **Correlation clustering** (`--clusterer cc`; the advisor's suggestion; Strainy itself uses karateclub's LabelPropagation). Read
  patterns are nodes weighted by read count. Edges follow the two-observation rule: patterns differing at < 2 shared columns are
  similar, at >= 2 dissimilar. The partition minimizes disagreeing pairs within clusters + agreeing pairs across. Solver: a
  deterministic pivot start, then local moves. No new constant.
  - P1 (harness ACC=a13 A=p23): D right **5,949** (81.5% of C) >= 5,842: **PASSES**; false moves 20 / 41,727 = 0.05%.
  - P2: **1 flag in 1 / 53** families (GWFAM163, 6 reads, c_unmatched): passes. Putative still finds GWFAM175's real copy
    (167 reads).
  - **Against label propagation (v6):** label propagation is slightly better on gorilla development data (P1 6,129 vs 5,949;
    P2 0 vs 1 false flag). Correlation clustering is better when the sample is another individual than the reference (the
    one-copy sweep: the erased copy's two alleles are no longer split, so it is detected down to 5 reads per allele instead of
    10). Both pass development. The choice is left to a held-out where both arms are registered.
- **Default = correlation clustering** (6e7ed568), chosen after the P5 held-out (both arms passed, correlation clustering better);
  `--clusterer lp` = v6 as registered.
- **Splice structure in the putative tier** (2026-10-23; user: "are we at any point using the intron/exon structure"). Until now the
  port used SNV columns only; introns were skipped.
  - **Diagnosis of P5's misses** (looked at after the verdict):
    - The 6 "no pile" deleted copies cannot be helped by structure: their reads never touch the survivor (most are unmapped,
      the rest align near the old locus).
    - The 2 free-slot copies and the 1 unclustered copy do splice differently from the survivor's own reads.
  - **First rule, rejected on development:** junctions the slot copy's other holders span but never use.
    - P1: +1 true, +2 false flags. P2: +1 false.
    - Structural contrast is symmetric: the survivor's allele looks as "foreign" as the deleted copy's.
  - **Rule kept:** a slot holder with >= MIN_READS of its majority junctions non-canonical on the slot copy's DNA (not GT-AG /
    GC-AG / AT-AC on either strand). Every isoform of a copy splices at its own sites, so alternative splicing cannot raise it.
    - P1 adds 2 true flags (GWFAM490 92 reads, GWFAM407 38) and 0 false.
    - P2 adds 0.
    - P5, post hoc (the threshold of 2 instead of >= half was set after GWFAM101's 7 of 17): catches GWFAM101's deleted copy
      (249 reads); 0 on the control.
  - It flags, it never certifies. `--no-structure` turns it off.
- **O2 v2** (after P7 failed): a tied read is assigned only when its haplotype was weighed for every copy the read ties to (the
  haplotype's options cover all the read's piles) and the edge is forced. In GWFAM328 (spent), a haplotype seen only on copy :10
  was trivially forced there, yet its tied reads also fit :11.
  - Development: 0 wrong everywhere.
  - Indistinguishable ties all abstain (GWFAM415, 407, 28, 401, 328).
  - The erased-copy simulations keep their fixes (163 / 163, 156 / 163, 43 / 48 placed on the certified new copy, 0 wrong).


---

## PREREG_o3_port_v6_heldout_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_port_v6_heldout_2026-10-23.md`)*

## Pre-registration: the O3 port v6 on a fresh gorilla fibroblast held-out (P5), 2026-10-23

**Status:** committed before the 20 families are drawn; no data of these families has been aligned or scored.

### Why
- **Development.** Port v6 passes both development panels (`docs/PREREG_o3_port_2026-10-23.md`, development log):
  - P1: D right 6,129 = 83.9% of C; false moves 0.02%.
  - P2: 0 of 53 families with a false flag.
- **Held-out so far.** P3 was spent on v4, and every change since was made on P1 / P2, so v6 has no held-out evidence.
- **Correlation clustering** passes development too (P1 5,949, P2 1 / 53) and is better in the one-copy simulation when the sample is
  another individual.
- **User constraints** (2026-10-23): "lets stay in gorilla, try the held out". Earlier: keep it small ("this might be overly
  complicating").

### Panel (fixed before drawing)
- **Families:** 20 drawn with `random.Random(1).sample` from P4's registered draw (`/mnt/linuxdisk/tmp/rna_allele/p4/panel.json`).
  P4's draw has 60 two-copy fibroblast families in neither the 53 development families nor refabsent's 34, and was never aligned
  or scored (`docs/PREREG_o3_two_tier_p4_2026-10-23.md`, "Not run").
- **Deleted copy:** P4's (Amendment 7's rule).
- **Reads:** P4's scored reads of the 20 (`link_test.py reads`: <= 500 primaries per copy, seed 1, the fibroblast BAM).
- **Deletion arm:** the `_pri` with only these 20 copies masked.
- **Control arm:** the unmasked `_pri` with every copy.
- **Alignment command:** Amendment 7 / 9's.

### Code (frozen)
- **Port:** `bench/o3_graph/o3port.py` as of e11cf8c6.
  - **Primary arm:** v6, the default options (label propagation).
  - **Secondary arm:** `--clusterer cc`.
  - Both use the stage's nets and products of the same arm (IsoCon side).
- **Stage:** `candidate_copies` frozen at `p4/bin/candidate_copies` (sha1 2e3cd592). It reproduces Amendment 15's products byte for
  byte on development batch g0. Batches of 10 families.
- **Scorers:**
  - Deletion arm: `merge_test.py score`, arm M.
  - Control arm: `control_test.py` lift / classify, Amendment 9's a / b / c / pri rule.

### Decision rules (the primary arm)
- **H5-1, sensitivity (deletion arm).** D right among the attainable D reads >= 0.50 x the attainable D reads.
  - Attainable D reads are truth-free, as in Amendment 13b: a record on a surviving copy in R.bam, or attributed into their own
    family's net.
  - Development P1: label propagation 4,994 / 8,663 = 0.58; correlation clustering 0.56.
- **H5-2, false moves (deletion arm).** <= 5% of the S reads. Development: 0.02%.
- **H5-3, specificity (control arm).** Families with >= 1 certified flag in class b, c or pri <= 3 of 20 (15%, P2's 8 / 53).
  Development: 0 / 53.
- **Verdict.** O3 v6 passes the held-out iff H5-1, H5-2 and H5-3 pass.
  - The secondary arm gets the same rules and is reported.
  - If only the secondary passes, that is a lead for a further held-out, not a switch made after looking.
  - The run is made once; nothing is changed after looking.

### Reported, not decided
- Look-harder flags by class (control) and by source (deletion).
- The landing loci.
- Each deleted copy's outcome: certified, look harder, in a slot, or no pile.
- Runtime.

### Result (run once, 2026-10-23; port at e11cf8c6, stage sha1 2e3cd592; work dir `/mnt/linuxdisk/tmp/rna_allele/p5`)
**Panel.** The draw is GWFAM101, 11, 122, 127, 148, 166, 196, 267, 279, 283, 304, 307, 357, 396, 408, 424, 491, 58, 75 and 80.
- 13,608 scored reads: 7,343 D, 6,265 S.
- Attainable D reads: 3,499 (2,381 on the survivor, 1,118 attributed).
- The stage's nets were replicated exactly on both arms (20 / 20 families).

| rule | primary: v6, label propagation | secondary: correlation clustering |
|---|---|---|
| H5-1 D right among attainable >= 0.50 | 2,308 / 3,499 = **0.660, PASS** | 2,901 / 3,499 = **0.829, PASS** |
| H5-2 false moves <= 5% (deletion arm) | 0 / 6,265, **PASS** | 0 / 6,265, **PASS** |
| H5-3 families with a false flag <= 3 / 20 (control) | 2 (4 flags: 2 pri, 1 b_allele, 1 c_unmatched), **PASS** | 1 (1 pri flag), **PASS** |
| D right, all D reads | 2,365 | 2,962 |
| deleted copies with >= 1 D read given back | 10 / 20 | 11 / 20 |

**Verdict: the O3 port v6 PASSES the held-out.** It is the first version to pass development (P1, P2) and a held-out.
- Correlation clustering passes as well and is better on every row here, as in the one-copy simulation.
- Development slightly favoured label propagation (P1 6,129 vs 5,949).
- Making correlation clustering the default would be a choice made after seeing P5. P5 then holds as evidence that each arm
  passes its own pre-registered rules, not for the choice between them.

**Scope.** These are gorilla KB3781 fibroblasts, two-copy families, and a deleted copy simulated by masking it in the reference.

**Each deleted copy's outcome** (reported, not decided; the D-majority haplotype's fate in the port's dump):

| outcome | label propagation | correlation clustering |
|---|---|---|
| certified (Strainy side) | 4 | 3 |
| recovered by the IsoCon side | 7 | 7 |
| in the survivor's free allele slot (Theorem A) | 2 | 3 |
| reads in the survivor's pile but no D haplotype | 1 | 1 |
| no pile | 6 | 6 |

The main loss is "no pile", not Theorem A. Those deleted copies' reads never land within 10 AS on the survivor: they align
elsewhere or not at all, and the shipped stage did not flag them either.

*Terminology (2026-10-23, after the run): the "look-harder" flags above are now called putative copies (user's request). Only
the name changed.*


---

## PREREG_o3_putative_clip_p10_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_putative_clip_p10_2026-10-23.md`)*

## Pre-registration: O3 putative copies with the soft-clip signal on a fresh gorilla fibroblast held-out (P10), 2026-10-23

**Status:** committed before the panel is drawn; none of its families has been aligned or scored.
**Goal milestone 2** (`ThesisVault/Topics/Goal_O2_O3.md`). The user said "work on 2" after P6 failed recall (register 1380).

### What changed since P6
The putative tier now has a third signal, `foreign_clipping` (`o3port.py` e3597fbc, default on): >= half of a slot holder's reads
with a near-best record on the slot copy are soft-clipped there by >= 50 bp that are not a poly-A / poly-T tail, and >= 6 of them
clip at one boundary (same side, same 20-bp bin). The deleted copy's reads stop matching the survivor where its UTR or an exon
diverges; a copy's own alleles align end to end.

**Development (P5 and P6, both spent; looked at after P6's verdict):**
- P6 deletion arm: free-slot recall 2 / 13 -> 5 / 13; 7 / 7 putative flags D-majority.
- P6 control arm: 1 flag (GWFAM405, as before; nothing new). P5 control: 0 haplotypes meet the rule.
- The thresholds (50 bp, half, 6 at one boundary, 70% tail) were set once before looking at P6's numbers and not tuned.

### Panel (fixed before drawing)
Every family eligible by P4's rule (`o3port_p4.py eligible_families`) that is in none of P4's, P5's or P6's panels: 58 two-copy
families, sorted by name. Everything else as in P6: Amendment 7's deleted copy, `link_test.py reads` (<= 500 primaries per copy,
seed 1), the deletion arm masks only these 58 copies, the control arm is the unmasked `_pri`, the stage `p4/bin/candidate_copies`
(sha1 2e3cd592) gives the nets and the IsoCon side. Work dir `/mnt/linuxdisk/tmp/rna_allele/p10`.

**Port:** `bench/o3_graph/o3port.py` as of e3597fbc, default options.

### Decision rules (as P6)
- **HP1, precision (deletion arm).** True putative flags / all putative flags >= 0.50 (true = D-majority, or class a by Amendment
  9's rule).
- **HP2, specificity (control arm).** Families with >= 1 putative flag in class b, c or pri <= 15% of the panel (8 of 58).
- **HP3, recall (deletion arm).** Of the deleted copies whose largest D-majority haplotype (>= 6 reads) holds a survivor's allele
  slot uncertified, the fraction flagged putative >= 0.25.
- **Verdict.** Pass iff HP1, HP2 and HP3. HP1 or HP3 with nothing to count: underpowered (no pass). Run once; nothing changed after
  looking.

### Reported, not decided
Each signal's share of the flags (tie, splicing, clipping); certified flags; runtime.

### Result (run once, 2026-10-23; port e3597fbc; work dir `/mnt/linuxdisk/tmp/rna_allele/p10`)
- Reads 38,675 (19,571 D); the stage replicated on both arms (58 / 58); certified flags 38 (deletion arm), 3 (control).

| rule | result | bar | verdict |
|---|---|---|---|
| HP1 precision (deletion arm) | 6 / 6 putative copies D-majority (GWFAM119 11 reads, 139 88, 171 6, 349 95, 363 313, 74 142) | >= 0.50 | PASS |
| HP2 control families with a false putative flag | 0 (no putative flag at all on the control) | <= 8 | PASS |
| HP3 recall over free-slot deleted copies | 5 / 8 = 0.625 | >= 0.25 | PASS |

**Verdict: the putative tier PASSES.**
- Reported: every one of the 6 flags carries the clip signature (>= 59 % of its reads clipped at one boundary), so the soft-clip
  signal is what fires. Misses: GWFAM380 (454 reads, all clipped when every record is counted, but not on the near-best records
  the port uses), GWFAM259 and GWFAM409 (no clipping: the deleted copy aligns end to end on the survivor).


---

## PREREG_o3_putative_heldout_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_putative_heldout_2026-10-23.md`)*

## Pre-registration: O3 putative copies on a fresh gorilla fibroblast held-out (P6), 2026-10-23

**Status:** committed before the panel is drawn; no data of its families has been aligned or scored.
**Goal milestone 2** (`ThesisVault/Topics/Goal_O2_O3.md`): "copies that fit a free allele slot are caught at a pre-registered
precision".

### What is tested
**Putative copies** (`o3port.py`, `putative_copies`) are the uncertified tier: haplotypes that hold an allele slot, so RNA cannot
prove them extra (Theorem A), but where either:
- **(tie)** >= half of the haplotype's reads tie (within 10 AS of their best) to another copy too, while < half of the reads of
  the slot copy's other holders do; or
- **(splicing)** >= 2 of its majority junctions are non-canonical on the slot copy's DNA.

A putative copy needs >= 6 reads and a co-holder on its slot copy.

**Development** (P1 / P2; P5 post hoc for the splicing threshold):
- P1:
  - 11 deleted copies end up in a survivor's free slot; 5 of them are flagged;
  - the flags are 5 D-majority, 2 on GWFAM175's real paternal-only copy, and 1 false (GWFAM415).
- P2: putative flags only on GWFAM175's real copy.

### Panel (fixed before drawing)
80 fresh two-copy fibroblast families, all outside the 53 development families, refabsent's 34, and every family ever aligned or
scored:
- **the 40 families of P4's draw that P5 did not use** (P4's draw: `docs/PREREG_o3_two_tier_p4_2026-10-23.md`, never scored);
- **40 more**, drawn with `random.Random(1).sample` from the eligible families P4 did not draw:
  - eligibility is P4's rule: `intervals.tsv`, every copy listed, >= 2 copies, >= 20 clean reads per copy, clean span <= 200 kb,
    G3;
  - the list is sorted by family.

The rest is as in P5:
- the deleted copy is Amendment 7's;
- the reads come from `link_test.py reads` on this panel (<= 500 primaries per copy, seed 1);
- the deletion arm has only these 80 copies masked; the control arm is the unmasked `_pri`;
- the stage is `p4/bin/candidate_copies` (sha1 2e3cd592), and the port's nets and IsoCon side come from it.

**Port:** `bench/o3_graph/o3port.py` as of 29d63d5b, default options (correlation clustering, putative copies with both signals).

### Decision rules (putative copies)
- **HP1, precision (deletion arm).** True putative flags / all putative flags >= 0.50.
  - A flag is **true** when most of its reads come from the deleted copy (D-majority), or when its consensus classifies as class a
    (a haplotype-only locus: a real reference-absent copy) by Amendment 9's rule.
  - Development P1: 7 / 8.
- **HP2, specificity (control arm).** Families with >= 1 putative flag in class b, c or pri <= 15% of the panel (12 of 80).
  - Development P2: 0.
- **HP3, recall (deletion arm).** Of the deleted copies whose largest D-majority haplotype (>= 6 reads) holds a survivor's
  allele slot uncertified, the fraction flagged putative >= 0.25.
  - Development P1: 5 / 11, with both signals. Two-copy families have one survivor, so only the splicing signal can fire.
- **Verdict.** The putative tier passes iff HP1, HP2 and HP3 pass. HP1 or HP3 with nothing to count (no flag, or no free-slot
  copy) does not pass: underpowered.
- The run is made once; nothing is changed after looking.

### Reported, not decided
- The certified tier on the same panel: certified flags, their classes on the control, and each deleted copy's outcome. No
  realignment scoring.
- Runtime.

### Result (run once, 2026-10-23; port 29d63d5b; work dir `/mnt/linuxdisk/tmp/rna_allele/p6`)
- **Reads and nets.** 49,778 reads (23,011 D). The stage's nets were replicated on both arms (80 / 80 families).
- **Certified flags** (reported): 40 on the deletion arm, 4 on the control arm.

| rule | result | bar | verdict |
|---|---|---|---|
| HP1 precision (deletion arm) | 2 / 2 putative copies are D-majority (GWFAM155 266 reads, GWFAM416 485) | >= 0.50 | PASS |
| HP2 control families with a false putative flag | 1 (GWFAM405, c_unmatched) | <= 12 | PASS |
| HP3 recall over free-slot deleted copies | 2 / 13 = 0.15 | >= 0.25 | **FAIL** |

**Verdict: the putative tier DOES NOT PASS** (not sensitive enough). What it flags is right (2 / 2) and it rarely flags wrongly
(1 / 80 families on the control). But in two-copy families only the splicing signal can fire, and 11 of 13 deleted copies that
sit in a survivor's free slot splice canonically there. Milestone 2 stays open.


---

## PREREG_o3_two_tier_p4_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o3_two_tier_p4_2026-10-23.md`)*

## Pre-registration: O3's two-tier output on a fresh gorilla fibroblast held-out (P4), 2026-10-23

**Status:** committed before the panel is drawn or built; no P4 data has been looked at.

### Why
- **The port does not pass alone.** Its best version (v5, `docs/PREREG_o3_port_2026-10-23.md`, development log) fails P1 by 80
  reads: D right 5,762 = 78.9% of C, bar 80%. Its precision is far better than the shipped stage's: false moves 0.02% vs 0.77%,
  false-flag families 0 / 53 vs 10 / 53.
- **The advisor's framing.** O3 "does not need to be definitive ... should suggest when to look harder". So O3 reports two tiers:
  - **Tier 1, certified:** the port v5's flags. These are Strainy-side certificates (more haplotypes than allele slots) and
    IsoCon-side candidates built mostly from reads in no phased pile.
  - **Tier 2, look harder:** the shipped stage's (`candidate_copies`) other flagged candidates.
- **Substrate.** Gorilla KB3781 fibroblasts only. The user, 2026-10-23: "we need to stick to fibroblasts in gorilla".

### Panel (fixed before drawing)
- **Source:** `intervals.tsv`, the source of Amendment 7's panel.
- **Families:** in neither the 53 development families nor refabsent's 34; every copy listed (rows = `n_copies`); >= 2
  copies; >= 20 clean reads per copy; clean span <= 200 kb per copy; G3 (the deleted copy overlaps no other gene).
- **Deleted copy:** the one last by (chrom, start), Amendment 7's rule.
- **Draw:** 60 families, `random.Random(1).sample` over the sorted eligible list (all of them if fewer).
- **Reads:** `link_test.py reads` (<= 500 primaries per copy, seed 1, the fibroblast BAM).
- **Arms:**
  - **P4d:** the masked `_pri`.
  - **P4c:** the unmasked `_pri`.
  - Both use Amendment 7 / 9's alignment command: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.

### Code (frozen)
- **Port:** `bench/o3_graph/o3port.py` as of 555e77db, default options (v5), `--isocon` = the stage's products of the same arm.
- **Stage:** `candidate_copies`, binary frozen at `p4/bin/candidate_copies` (sha1 2e3cd592). It reproduces Amendment 15's stage
  products byte for byte on development batch g0. Default options, batches of 10 families.
- **Scorers:**
  - Deletion arm: `merge_test.py score`, arm M.
  - Control arm: `control_test.py` lift / classify (Amendment 9's a / b / c / pri rule) / score (arm C).

### Decision rules (tier 1)
- **H4-1, specificity (P4c).** Families with >= 1 tier-1 flag in class b, c or pri <= 9 of 60 (15%).
  - Development: P2 0 / 53.
  - v4 on P3: 3 / 30.
- **H4-2, false moves.**
  - P4d arm M: <= 5% of the S reads.
  - P4c arm C: <= 5% of all scored reads.
- **H4-3, sensitivity (P4d).** Tier-1 D right among the attainable D reads >= 0.50 x the attainable D reads.
  - Attainable D reads are truth-free, as in Amendment 13b: a record on a surviving copy in R.bam, or attributed into their own
    family's net.
  - Development P1: 4,627 / 8,663 = 0.53.
- **Verdict:** the certified tier holds on fresh fibroblast data iff H4-1, H4-2 and H4-3 all pass. The run is made once; nothing
  is changed after looking.

### Reported, not decided
- **Tier 2** = the stage's flagged candidates that are not in tier 1 (sharing < half of their members with every tier-1 flag):
  - its false-flag families on P4c;
  - the stage arm's own D right among the attainable reads, and its false moves on P4d.
- **Per deleted copy:** caught by tier 1, by tier 2 only, or by neither.
- **Also:** a per-family table, and the runtime.

### Not run (2026-10-23, the same day)
The user stopped the build before any P4 data was scored: "I think this might be overly complicating, can we somehow only focus on
one erased copy? We could even simulate it so it's clean". Only the panel draw, the masked genome and the scored reads were made
(`/mnt/linuxdisk/tmp/rna_allele/p4`). Nothing was aligned or scored. The work moved to `docs/O3_ONE_COPY_SIM_2026-10-23.md`.

---

## Copy phasing in Rust (`copy_phasing`, the port of `o3port.py`), 2026-10-24

Plan: `docs/superpowers/plans/2026-10-24-o3port-rust.md`. The Rust binary is named for what it does (user, 2026-10-24): it phases each
copy's reads into haplotypes, matches haplotypes to allele slots, certifies or flags copies the reference lacks, and assigns tied reads.
Product file names (`.candidates.tsv` ... `.o2.tsv`) are unchanged, so the harness and scorers read either port.

**Byte-identity** (`bench/o3_graph/o3port_golden.sh all`, cases in `o3port_golden.tsv`): 78 of 78 cases identical in every product, the
`--dump` JSONL and the per-family progress lines, at `--threads 1` and at `--threads 4`. Cases: 9 o2p14 dirs × 7 flag sets (default,
--semi, +--guard, +--no-reconstruct, +--anchor, +--one-site, +--prefer), every other switch alone on one dir, o2p12 GWFAM176, P5 and P10
both arms with the IsoCon side, the refabsent `h` panel (--semi --guard; --isocon-rule via-b). The Python is itself deterministic: the 72
small cases are identical under PYTHONHASHSEED 1 and 2. Numerics that had to follow CPython: `math.lgamma` (CPython's own Lanczos, not
glibc: they differ on half the inputs), `sum()` of floats (Neumaier since 3.12), `repr(float)` in the dump, `'%.3f'`.

**Speed** (P10 control arm, `R0.bam`, 40 families, IsoCon side): Python 166.4 s / 1.9 GB; Rust 19.1 s / 0.69 GB at 1 thread, 5.4 s /
1.16 GB at 4 threads.

**The nets step in the stage** (`candidate_copies --net-reads PATH`): the per-read net table `accept_o3_candidates.py nets` writes
(`family read used role via class in_reads_tsv`, `role` NA: truth labels stay in the harness). Not in the stage key (it adds a product);
stored in the cache entry and replayed only when asked for. Against the Python nets' `net_reads.tsv` (every column but `role`), one batch
each: P5 g0 8,166 rows, P10 g3 6,632, refabsent `h` g2 4,572 (19 pass-B poorly placed reads): identical; the batch's 9 stage products equal
the frozen binary's. **Switch:** `PORT=rs bench/o3_graph/o3port_p4.sh port <arm>` runs `copy_phasing` into `p<arm>.rs` (on the nets step's
table when present, else on the stage's `--net-reads` tables, `NET_READS=1` at `stage`); on P10 `port c` its 10 products equal the
frozen Python `pc` products, 8.8 s. Register row 1402.

**Final review fixes (2026-10-24).** (1) A haplotype's copy options were in `--copies` order, not sorted as in the Python. `b_matching` takes
the first free option, so the slot copy could differ whenever a copies file lists a copy before a lower-sorting one (P5 / P10 list `:1`
before `:0`). The 78 golden cases passed only because their one two-option haplotype was re-sorted on another path. Fixed; 3 golden cases
with the copies reversed (`p14_d1_control_rev_*`) added. **81/81 identical at 1 and 4 threads.** (2) The `PORT=rs` concat branch had a
Python syntax error. Fixed; run on a P10 stand-in (per-batch tables, no nets-step file): the table and all 10 products equal the frozen
`pc` products. (3) The `candidate_copies` stage key moved to v3, so a pre-`--net-reads` cache entry is never reused.
