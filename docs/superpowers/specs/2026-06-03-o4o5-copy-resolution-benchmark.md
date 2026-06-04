# O4/O5 copy-resolution benchmark — the discriminating multi-copy benchmark

**Date:** 2026-06-03. **Status:** built + measured. **Harness:**
`bench/tandem_attribution/o4o5_copy_benchmark.sh` (multi-seed). **Scorer:**
`bench/tandem_attribution/score_attribution.py`. **Generator:**
`bench/tandem_attribution/gen_synthetic.py` (extended: `--distinct-isoforms`, `--starve-copy`).

## Why this benchmark exists

The audit (2026-06-03) found the existing synthetic O4 fixture non-discriminating: copies are
spatially separable, so the **non-VG baseline scores the same 100/100 as VG** on gffcompare
intron-chain Sn/Pr. The benchmark was supposed to be a fixture where VG beats the baseline.

A grounded discrimination study (workflow `o4o5-benchmark-discrimination`, 3 adversarial
investigations) established a **harder, more honest truth**:

> **gffcompare intron-chain Sn/Pr is the WRONG metric — VG cannot beat the baseline on it.**
> When every copy has full-length reads the baseline scores 100/100 (one transcript per
> separated bundle) and VG can only TIE or, via secondary-alignment intake, LOSE precision.
> Verified across copies 2–6, spacings 13k–30k, with/without starvation, with/without
> `RUSTLE_VG_TANDEM`. No synthetic regime makes VG win on chain Sn/Pr. The tandem decomposer
> only ever fragments correct calls into false positives (`vgt_copies_base_lacks=[]` in 6/6
> fired cases); a starved copy whose 5′-most junction is observed by no read is unrecoverable
> by baseline, VG, AND tandem alike.

The real discriminator is **copy attribution**: which paralog copy did each ambiguous
(multi-mapping) read come from? VG's fingerprint-EM answers this against per-read ground truth
(reads named `c{src}_r{j}`). The baseline has **no copy concept at all** (zero
`copy_id`/`copy_confidence`/`family_*` attributes), so the metric is **UNDEFINED** for it — not
merely worse. This benchmark demonstrates a *capability* VG provides and the baseline does not,
and characterizes the identifiability boundary including the critical honesty check: does VG
**abstain** when copies are identical (the DAZ non-identifiable limit), or fabricate?

## Fixture

A merged-but-separable tandem pair: `--copies 2 --distinct-isoforms --spacing 16000` (>
gene_len 11600, so copies do not physically overlap) `--identity <id> --reads-per-copy 80`.
Every read multimaps to both copies (`minimap2 -ax splice:hq -uf -N20 -p0.5`), so the ambiguous
reads are genuine. `--distinct-isoforms` gives each copy a copy-specific cassette-exon skip
(distinct intron chains). reads-per-copy 80 makes the family clear the **default**
`--vg-family-min-shared=10` so Part 1 runs at default settings.

## Results (measured, 5 seeds {1,3,5,7,11})

### Part 1 — copy attribution, DEFAULT settings (identity 0.97)
| metric | VG | baseline |
|---|---|---|
| copy-attribution accuracy (multimappers) | **1.00** (5/5 seeds) | **undefined** (0 copy attrs) |
| decisive accuracy | **1.00** | undefined |
| chain Sn/Pr | 100/100 (4 seeds), 50/50 (1 seed) | 100/100 |

VG attributes every ambiguous multimapper to its true source copy; the baseline produces no
copy metric. Chain Sn/Pr confirms it does **not** discriminate (baseline ties or beats VG; the
one 50/50 is the family-dropped duplicate-FP artifact).

### Part 2 — identifiability spectrum (relaxed `--vg-family-min-shared 2`, documented)
Relaxed because the real-data spurious-family guard (min_shared=10) drops this small KNOWN
family at high identity; relaxing lets the EM run so its behavior across the spectrum is
measurable. The family is controlled and known-real, so the guard is not needed here.

| identity (exonic SNPs) | acc | dec_acc | dec_frac | abstain | reading |
|---|---|---|---|---|---|
| 1.0 (0 SNP) | 0.58 | 0.00 | 0.00 | **1.00** | **abstains 100% — fabricates nothing (DAZ limit) ✓** |
| 0.999 (1 SNP) | 0.55 | **0.44** | **0.75** | 0.25 | **OVERCONFIDENT — 75% decisive at chance-level accuracy ✗** |
| 0.99 (≈9 SNP) | 1.00 | 1.00 | 1.00 | 0.00 | identifiable → correct |
| 0.97 (≈27 SNP) | 1.00 | 1.00 | 1.00 | 0.00 | identifiable → correct |
| 0.95 (≈45 SNP) | 1.00 | 1.00 | 1.00 | 0.00 | identifiable (1/5 seeds: family dropped) |
| 0.90 | — | — | — | — | below ambiguity: reads map uniquely, no family forms (5/5 dropped) |

## Honest findings

1. **The metric is copy attribution, not chain Sn/Pr.** On the merged fixture VG provides
   correct per-read copy resolution (100% in the identifiable band); the baseline structurally
   cannot (no copy concept). This is the load-bearing O4/O5 result the old 100/100 fixture
   lacked — but it is a *capability demonstration*, not a head-to-head Sn/Pr win.
2. **VG abstains at the non-identifiable limit (id 1.0): 0 decisive calls, fabricates
   nothing.** This is the DAZ3-false-positive discipline, demonstrated on a controlled fixture.
3. **VG is OVERCONFIDENT at the 1-SNP boundary (id 0.999): `dec_acc 0.44` at `dec_frac 0.75`.**
   The abstention boundary is too sharp — at exactly-identical it abstains, but at 1 diagnostic
   SNP it makes confident calls that are worse than chance. This is the benchmark's main
   actionable finding: an EM-calibration gap at the boundary, the exact regime where
   fabrication risk lives. **Surfaced, not hidden.**

## Limitations (honest)

- **Synthetic**, truth = read-name prefixes (same caveat as the O3/O4 oracle). Does not
  establish real-data (GGO) copy-attribution accuracy — that remains the cross-cutting gap.
- Part 2 relaxes the family filter (documented). Part 1 runs at default settings.
- Discovery is stable in id 0.95–0.99 with rpc 80; at id≥0.999 it fluctuates near the filter
  threshold (a small-fixture artifact, not a model property).
- The benchmark scores **attribution**, not per-copy abundance; abundance calibration on a
  genuinely two-expressed-copy family is a separate follow-up.

## Next steps

- Wire Part 1 (id 0.97 copy-attr=1.0) + the id 1.0 abstention into `run_oracle.py --check` as
  regression guards.
- Fix the id 0.999 overconfidence (EM should abstain, or widen uncertainty, at 1 diagnostic
  site) — then re-measure.
- The real-data analogue (RBMY / a genuinely two-expressed paralog pair with curated per-read
  truth) is the validation the audit named; this synthetic is the controlled precursor.
