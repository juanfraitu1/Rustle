# Pre-registration: O2's regimen on almost-tied reads, in simulation and on real reads (round 24 + A119b fragments, 2026-10-19)

**User (2026-10-19, /goal):** "establish a working O2 regimen from both simulations and real reads, O2 must get almost tied
multi-mapping reads and help assign via PSVs from which copy did they actually come from".

## The regimen (`bench/o3_graph/o2psv.py svp`, code 26262631, defaults: k = 3, M = 10, fit 0.1%)

- **Domain (almost tied):** a read whose alignments within 10 AS of its best (Eichler) lie on two or more copies of one
  family and on no locus off it. Each alignment belongs to one copy: its strand, then exon overlap.
- **Evidence (PSVs):** per tied copy, the substitutions of its best tied alignment. Clips and indels are left out.
- **The read is placed** when either:
  - **two PSVs:** one copy has at least 2 fewer substitutions than every other; or
  - **one sample-validated PSV:** the read carries copy X's base at a PSV column, and the sample's own reads placed on the
    rival copy Y (AS margin >= M = 10 over every other copy: Eichler's well placed; the read itself left out) cover Y's
    position at least k = 3 times
    and none shows that base. The one copy no validated PSV excludes takes the read, provided the read's bases contradict
    it at no PSV column.
- **Otherwise** the read keeps the set of copies neither rule removes; disagreeing rules abstain. A read its best copy does
  not explain is O3's (`psv_unexplained`): more than 0.1% substitutions, HiFi's Q30 rate.

Why both rules:
- in one genome, almost-tied reads differ by at most about one PSV, which is enough there (round 23: 2,045 right, 0 wrong);
- on a real donor, a single PSV is the donor's own allele often enough to matter (HG002: one-PSV placements 93.4% right).

The sample-validated PSV is the bridge: it is used only where the sample's own copy shows it holds.

## Development (spent)

**HG002, donor truth:** all reads within 10 AS on >= 2 copies, poor_match unscored.

| HG002 setting | right | wrong | precision |
|---|---|---|---|
| **the regimen** | **3,169** | **35** | **98.9%** (37 own_tie placements) |
| two PSVs alone | 2,635 | 25 | 99.06% |
| validated alone | 534 | 10 | 98.1% |
| fit at the sample's own p90 / p99 instead of 0.1% | | | 98.3% / 97.9% (more placed) |
| validator margin 1 instead of 10 | | | the same |
| minimap2's primary on the same reads | 7,091 | 3,097 | 69.6% |

**Simulation:** three spent CHM13-only round-23 worlds.
- 124 placed, 0 wrong; every set holds the true copy (1,227 of 1,227); copies 2.01 -> 1.85-2.00.
- minimap2's primary is 55-64% right on the same reads.
- The fixed 0.1% fit costs power under the noisier public profile: 178 reads unexplained in one world.

## Held-out test A: simulation, one genome

- **Families:** round 24's draw: 21 fresh families, 90 copies (`families24.json`, set `r24`, sha1 8611cffce41b; the pool
  less round 23's).
- **Worlds:** INTACT, CROSS and SILENT, under `profiles/{A119b,public}_testis_v4_ref.json` (no donor alleles), seeds 1401,
  1402 and 1403. That makes 6 cells and 18 worlds.
- **Run:** `fastset.py run` with round 22's `$V17` builds the worlds. Then, per world, `o2psv.py copies-world` and
  `o2psv.py svp E0.bam E0.bam copies genome.ref.fa --k 3` (the defaults), scored against
  `reads.truth.E0.tsv`.
- **Bars (pooled):**
  - **A1:** wrong at most 0.5% of placements; floor 200 placements.
  - **A2:** sets hold the true copy for at least 99% of the domain's reads.
  - **A3:** O2's precision is above minimap2's primary's on the same reads.

## Held-out test B: real reads, A119b testis (a donor the regimen was not developed on)

- **Sources:** per pool family (`r19seq`, 439 families), 60 reads by a seeded hash (`o2frag.py sources`, seed 1). Each is
  a MAPQ-60 primary on a copy (strand and exon overlap), at least 1 kb long. Scripts:
  `fast/heldout24/run_sim.sh` and `o2frag_a119b/run_real.sh`.
- **Fragments:** 3 per source, 300-1,200 bp, seeded (`o2frag.py fragments`). Sources and fragments are realigned to whole
  CHM13 with A119b's own parameters: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.
- **Truth:** a fragment's truth is its source's copy, when the realigned full read lies on that copy at least 30 AS clear of
  any other locus (`o2frag.py truth`). On HG002, reads that clear agree with the donor's assembly 98.4% of the time, so
  this proxy carries 1-2% noise.
- **O2:** `o2psv.py svp FRAG.bam A119b.t2t.bam pool_copies target.fa --names truth.tsv --exclude sources.tsv` (a source
  never validates its own fragment), same k, M and quantile.
- **Bars (pooled):**
  - **B1:** wrong at most 2% of placements; floor 100.
  - **B2:** O2's precision is above minimap2's primary's on the same fragments.
  - **B3:** sets hold the true copy for at least 95%.

## Verdict

- **PASS:** A and B pass, each powered.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

## Not claimed

- That identical copies can be split: they stay a set.
- Copies the reference lacks (O3).
- Truth beyond B's proxy. HG002's donor truth was development.

**Rules:** code frozen at 26262631. A bug is recorded as an amendment before any rerun.
