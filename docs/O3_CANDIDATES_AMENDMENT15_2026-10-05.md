# O3 candidates, Amendment 15/15b — the consensus-defect fix and its re-runs, 2026-10-05

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendments 15/15b (written 2026-10-03, before any re-run).
Defect diagnosis: `docs/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md` (byte-identical reproduction, 33 clusters of 5 families).
Code: commit `a13b817f` (the correction + regression tests); harness commit `209a5fd7` (`ACC=a15h`, env-overridable work dirs).
This doc records the implementation and every registered re-run; nothing here re-tunes a registered rule.

## The correction (as registered, no new constants)

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

## Re-runs and verdicts (same bars, nothing re-tuned)

All runs with the corrected binary (`a13b817f`, sha1 6fc286a0), the registered substrates, no `RUSTLE_CACHE_DIR` (every
minimap2 call computed). Work dirs `a13a15/`, `a14a15/`, `a15h/` — the registered `a13/`, `a14/` products untouched.

### A13 again (the 53-family deletion held-out) — A13-1/2/3 all PASS

- Stage: 53 families, 363 clusters, 77 candidates, **49 flagged** (vs 82 flagged in the uncorrected run), 1,802.9 s = **30.0 min** (bar ≤ 40).
- **A13-1 PASSES**: D right **7,882** ≥ 0.80 × C = 5,842.4 (107.9 % of C; identical to the uncorrected run's 7,882 — the
  correction lost no correct placements); false moves 322/41,727 = **0.77 %** ≤ 5 %.
- **A13-2 PASSES**: the unions keep **99.33 %** (pooled) / 99.43 % (isolated) of the components' reads (bar ≥ 95 %).
- The flags' composition changes exactly as the diagnosis said it should: D-derived **30/30 kept**, survivor-derived
  **46 → 13** (the 13 kept carry the divergent-but-genuine flags: their best-hit identity×coverage was already high in
  the uncorrected run), "elsewhere" 6 → 6. Deleted copies with ≥ 1 D-derived flag: 25/53 (unchanged); the cause table is
  unchanged in shape (23 no-D-reads-in-net / 2 no-cluster / 3 linked).

### A14 again (the no-deletion control) — C2' HOLDS, **C1' FAILS**

- Stage: 53 families, 415 clusters, 40 candidates, **17 flagged in 10 families** (uncorrected: 54 class-c flags).
- Classes: 1 a_haplotype_only (GWFAM175), 2 b_allele, 14 c_unmatched.
- **C1' FAILS**: families with ≥ 1 class b/c/pri flag = **10/53 = 18.9 %** (bar ≤ 8 = 15.1 %). The 16 false flags are the
  defect doc's *other* category, not the fixed mechanism: median 11.5 reads, mismatch-dominated best hits
  (median whole-length d 0.036, median id×cov 0.974) — few-read or paralog-mixed consensus, conservative by construction.
- **C2' HOLDS**: false moves 317/59,013 = **0.54 %** ≤ 5 % (295 out of 317 from GWFAM175 alone).

### Held-out H (Amendment 10's read set, nothing deleted; the 30 disjoint families) — H3 HOLDS, **H1 FAILS**

New `ACC=a15h` harness mode (commit 209a5fd7): `refabsent/R0.bam` (32,219 reads), W.copies restricted to the 30 disjoint
families (66 copies; GWFAM175 is dev-overlapping and correctly absent), classification by Amendment 9's rule.
- Stage: 30 families, 208 clusters, 18 candidates, **11 flagged in 6 families**; 507.7 s = 8.5 min.
- Classes: 1 a_haplotype_only, 2 b_allele, 8 c_unmatched (identical at δ/2, δ, 2δ).
- **H1 FAILS**: 6/30 families with ≥ 1 class b/c/pri flag (bar ≤ 4/30).
- **H3 HOLDS**: false moves 73/27,757 = **0.26 %** ≤ 5 %.

### Decision rule (Amendment 15, registered): the default-on flip requires A13-1/2/3 AND C1'/C2' AND H1/H3

C1' and H1 fail → **the stage STAYS OPT-IN** (`--candidates`, ruling R14 stands); the corrected stage replaces the
uncorrected one as the opt-in arm — the duplication defect is gone (below), and the remaining false flags are a named,
different phenomenon (few-read consensus) for a future amendment, not this one.

## Reported beside (as registered)

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

## Files

| where | what |
|---|---|
| `a13a15/` | the corrected A13 run (all steps, incl. arm M, comparator, keep, report) |
| `a14a15/` | the corrected A14 run (classify, arm C, creport) + `wholebam/` R23 probe |
| `a15h/` | the held-out H run (classify, arm C) |
| `/tmp/dupsig.*.names` | the per-union off-diagonal sets behind the signature counts (session scratch) |

## Corrections (2026-10-06)

Found by an adversarial re-check of the Locus Anatomy page against the run products; each recomputed here from the products. The verdicts above (A13-1/2/3 PASS, C1' and H1 FAIL) do not change.

1. **A13-1, D right.** The uncorrected run's value is **7,898** (108.1 % of C) with 323/41,727 false moves (`a13/comparator.out`, `a13/score.out`), not 7,882. The correction cost 16 correct placements (-0.2 %) and one false move. "Identical to the uncorrected run's 7,882 — the correction lost no correct placements" (A13 section, first bullet) is wrong, and so is "unchanged by the fix" in register row 1247.
2. **A13-3, time.** The registered measure is the sum of the five batches' `/usr/bin/time` Elapsed: corrected **2,002.1 s = 33.4 min** (205.88 + 330.65 + 513.44 + 557.72 + 394.39), uncorrected **1,555.8 s = 25.9 min**: the majority test costs **+29 %**, not "~15 %". The 30.0 min above is the stage's own monotonic clock (1,802.9 s), which runs about 11 % behind Elapsed. The bar (40 min) is met on either clock.
3. **GWFAM175_0 and GWFAM331_0 are not few-read flags.** GWFAM175_0 has 183 reads (it is the class-a true copy, GWFAM175_B0) in both the A13 and A14 corrected runs; GWFAM331_0 has 143 reads in A13 and 129 in A14 (`cand.candidates.tsv`). The "few-read flags counted in C1'" sentence of the duplication-signature bullet holds for neither.
