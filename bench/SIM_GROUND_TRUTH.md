# Fully-simulated ground-truth benchmark — the airtight, non-circular validation (2026-06-29)

Every prior accuracy number on real GGO data carries some circularity: the "silver" label is minimap2's own
primary, the copies are RNA-assembled, the PSV columns are RNA-discovered. This benchmark removes *all* of it.
We **plant** a 2-chromosome genome — positions, divergence, copy number, exon/intron structure — label every
read with its TRUE family/copy, run the unmodified pipeline, and check whether it recovers exactly what we
planted. Nothing is borrowed from the pipeline's own output, so there is no circular reference.

- `bench/sim_genome.py` — plants the genome + labelled reads (deterministic, seed 20260629).
- `bench/sim_run.sh` — end-to-end: plant → minimap2 → `gw_family_catalog` (O1) → `copy_assign` (O2) → score.
- `bench/sim_eval.py` — scores the pipeline output against the planted truth.

## What is planted (`simA` 198 kb, `simB` 198 kb, 920 reads)

| family | copies | regime it probes |
|---|---|---|
| **K0tandem** | 3 IDENTICAL tandem | the K=0 floor — every read must be certified TIED, never guessed |
| **ladder** | 4 tandem @ 0 / 0.3 / 0.8 / 1.5% | the resolvable frontier as divergence rises |
| **collapse** | 5 near-identical, ≈6 kb, only 6 PSVs | **collapsed segdup** — multimaps (MAPQ-0) yet PSV-resolvable: the regime the gate is *for* |
| **cnv** | 3 tandem, unequal expression (80/40/20) | abundance / quantification |
| **xchrom** | 2 on DIFFERENT chromosomes, 0.3% | cross-chrom family detection + assignment |
| single ×4 | single-copy genes | must NOT form a family (FP control) |
| domshare ×2 | two genes sharing ONE 150 bp exon | must NOT merge (over-merge FP control) |

## Results (all non-circular — the read name carries the truth)

### O2 — per-read copy assignment (`copy_assign`, the significance gate)

Scored over **multimapping reads** (minimap2 keeps a secondary / MAPQ-0 — the reads that actually need the
gate), conditioned on the true copy being present in the model (a copy the detector did not recover cannot be
assigned to — that is an O1 miss, not an O2 error, and is counted separately):

```
ACCURACY on assigned reads: 460/460 = 100.0%   (= the TRUE planted copy)
  div=0.000 (copy0 of a divergent family): 120/120 = 100%
  div=0.001 (COLLAPSED segdup, 6 PSV/6kb):  200/200 = 100%   <- the gate's regime, minimap2 cannot
  div=0.003:                                  40/40  = 100%
  div=0.005:                                  60/60  = 100%
  div=0.008:                                  40/40  = 100%
K=0 floor (3 identical copies): 120/120 = 100% TIED   (certified unresolvable — never a fabricated pick)
misassignments: 0
```

- **The collapsed-segdup row is the headline.** 5 copies that are 99.9% identical (6 PSVs spread over a 6 kb
  transcript) multimap with low/zero MAPQ — minimap2 cannot separate them — yet the gate, which scores *only*
  the PSV columns, assigns **all 200 reads to the correct copy**. This is the in-genome demonstration of the
  thesis claim, on planted ground truth.
- **The K=0 floor is exact:** the 3 byte-identical copies yield 120/120 TIED with `min_p = 1.0` — the
  identifiability certificate fires, so the gate abstains rather than guessing. Irreducibility is respected.

### O1 — family detection (`gw_family_catalog`, cross-chrom-aware)

```
K0tandem (3): RECOVERED 3/3      collapse (5): RECOVERED 5/5
ladder   (4): RECOVERED 3/4      cnv      (3): RECOVERED 3/3
xchrom   (2): RECOVERED 2/2  (cross-chromosome)
over-merged families (span >1 planted family): 0
single-copy / domain-sharer reads wrongly assigned to a family: 0
```

- **0 false merges:** the four single-copy genes and the two domain-sharers (sharing a 150 bp exon, <13% of the
  transcript) are correctly left out of every family — the conflict-graph definition does not over-merge.
- **Cross-chrom works:** the 0.3% cross-chromosome pair forms a single 2-copy family (`cross_chrom=true`).
- **The one non-recovery is the definition working as designed.** ladder copy3 (1.5% diverged) is *not* merged
  into the ladder family: at 1.5% its reads are no longer *confusable* with the other copies (no de-tie
  conflict edge), so it is correctly reported as a separate locus rather than a 4th ladder copy. The
  conflict-graph family is the set of mutually-confusable copies — and a copy that has diverged past
  confusability genuinely is not one.

## The conceptual finding this makes airtight

Across the whole genome, **the divergence that makes copies PSV-resolvable is essentially the same divergence
that makes them uniquely mappable (MAPQ>0).** Uniform divergence ≥~1% lifts minimap2's MAPQ the moment it is
large enough to create PSVs (the ladder 1.5% copy, an early-divergence cross-chrom pair). So the genuinely
*ambiguous* (multimapping) mass concentrates at two places:

1. **identical copies (K=0)** — where the gate correctly certifies TIED (the information limit), and
2. **collapsed segdups** — many near-identical *long* copies with a *few* concentrated PSVs, where a handful of
   mismatches do not move MAPQ but the gate still resolves every read.

This is exactly the regime split the read-level `sim5x` ladder and the DNA-supervised decode operate in, and it
explains the real-GGO result (the unassignable mass is the K=0 identical floor, not gate failure). The gate's
value is precisely the collapsed-segdup band that whole-read mapping cannot touch — and here, on planted truth,
it resolves that band at 100% while never guessing on the irreducible floor.

## Reproduce

```bash
bash bench/sim_run.sh        # plant → map → catalog → assign → score (deterministic)
```
Outputs land in `winloci_scratch/simgw*` ; the planted truth is `simgw_truth.tsv`.
