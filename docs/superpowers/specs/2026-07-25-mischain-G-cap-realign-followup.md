# Mis-chain rescue via targeted `-G`-capped re-alignment — follow-up

**Date:** 2026-07-25
**Status:** promising lead, NOT yet validated (do not merge). Parked for follow-up after advisor meeting.

## The finding

Compound-failure members like the NCF1 family (NCF1 / NCF1B / NCF1C) miss on the real BAM because their
reads **mis-chain**: minimap2 chains a read's two ends across a spurious giant intron to a near-identical
paralog copy, instead of keeping the read local. The local pool then falls below the seeding floor and no
consensus forms.

Evidence (chr7 NCF1 cluster, `soto_regions.bam`):

- Local coverage exists: 44 / 15 / 17 primary reads overlap NCF1 / NCF1B / NCF1C.
- The mis-chain "introns" are **234–387 kb** — larger than NCF1's real introns (<5 kb) *and* larger than
  minimap2's default splice cap (`-G 200000`).
- Re-aligning the 76 cluster reads to CHM13 with a max-intron cap:

  | `-G` | reads staying local | reads still mis-chained (>50 kb intron) |
  |------|--------------------:|----------------------------------------:|
  | 200000 (≈ current) | 34 / 76 | **42 / 76** |
  | 50000 | **76 / 76** | **0 / 76** |

  A 50 kb cap hands back **every** mis-chained read as a local alignment. At real coverage that restores
  NCF1/B/C's full 44/15/17 local reads — enough to seed **without** the topup simulation. This would be a
  genuine recall gain, not a "recoverable-with-ideal-coverage" asterisk.

Topup already corroborates the mechanism: it "recovered" NCF1 only by resampling its own reads to ideal
depth, i.e. by brute-forcing enough locally-chaining reads. The `-G` cap achieves the same local pool from
the *real* reads.

## Why this differs from the REJECTED mis-chain read-salvage (branch `mischain-salvage`)

The rejected approach **split** each mis-chained read into two fragments (local + paralog). The paralog
fragment created a phantom node that over-connected the conflict graph and collapsed neighboring families
genome-wide (net +0, 6 gained / 6 lost). A `-G`-capped re-alignment instead keeps the read as **one** local
alignment (the paralog end soft-clips off) — no second node, no graph over-connection. It attacks the
artifact at the aligner, not downstream.

## Why it is NOT a one-line default

`-G 50000` is aggressive: many real human genes have introns >50 kb (some >200 kb). A global cap would
fragment those genes' alignments. So the fix must be a **targeted re-alignment pass**: only reads carrying an
implausibly long intron *whose far segment lands on a paralogous copy of the same locus* get re-aligned under
the cap. Length alone is not the discriminator — "spans to a homologous copy" is.

## Two-step validation plan (before believing it)

1. **Focused seeding test.** Re-align the NCF1 cluster reads genome-wide with `-G 50000`, rebuild the mini-BAM,
   run detection on the NCF1 family bed, confirm NCF1/B/C seed at real coverage. (Proof of concept.)
2. **Genome-wide net check.** Same discipline that killed the read-split: run the full `--cross-chrom` Soto
   catalog with the targeted re-align pass ON vs OFF; require net-positive recall with no family regressions
   elsewhere. Report gains AND losses honestly.

## Fallback if the aligner fix is too blunt genome-wide

The user's **coverage-presence** idea: for a defined family-member locus, count spliced/unspliced reads; if
present, call the member detected. Mis-chain-immune, but (a) needs the member locus from a non-circular source
(clean family projection, not Soto's bed), and (b) detects *presence* only — the MAPQ-0 multimapper ambiguity
still blocks resolving *which* copy (O2 stays K=0). Valid as an O1 detection layer.

## Reproduction

- Local coverage / intron sizes: `samtools view -F2308 soto_regions.bam <locus>` + CIGAR `N`-scan.
- Re-align test: extract cluster reads (`samtools fastq -F2308`), `minimap2 -ax splice:hq -uf -N 50 -p 0.1 -G <G>`,
  count primary alignments with a >50 kb `N` in CIGAR. See `/home/juanfra/winloci_scratch/ncf1_realign/`.
