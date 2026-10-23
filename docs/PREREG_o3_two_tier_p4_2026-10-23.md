# Pre-registration: O3's two-tier output on a fresh gorilla fibroblast held-out (P4), 2026-10-23

**Status:** committed before the panel is drawn or built; no P4 data has been looked at.

## Why
- **The port does not pass alone.** Its best version (v5, `docs/PREREG_o3_port_2026-10-23.md`, development log) fails P1 by 80
  reads: D right 5,762 = 78.9% of C, bar 80%. Its precision is far better than the shipped stage's: false moves 0.02% vs 0.77%,
  false-flag families 0 / 53 vs 10 / 53.
- **The advisor's framing.** O3 "does not need to be definitive ... should suggest when to look harder". So O3 reports two tiers:
  - **Tier 1, certified:** the port v5's flags. These are Strainy-side certificates (more haplotypes than allele slots) and
    IsoCon-side candidates built mostly from reads in no phased pile.
  - **Tier 2, look harder:** the shipped stage's (`candidate_copies`) other flagged candidates.
- **Substrate.** Gorilla KB3781 fibroblasts only. The user, 2026-10-23: "we need to stick to fibroblasts in gorilla".

## Panel (fixed before drawing)
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

## Code (frozen)
- **Port:** `bench/o3_graph/o3port.py` as of 555e77db, default options (v5), `--isocon` = the stage's products of the same arm.
- **Stage:** `candidate_copies`, binary frozen at `p4/bin/candidate_copies` (sha1 2e3cd592). It reproduces Amendment 15's stage
  products byte for byte on development batch g0. Default options, batches of 10 families.
- **Scorers:**
  - Deletion arm: `merge_test.py score`, arm M.
  - Control arm: `control_test.py` lift / classify (Amendment 9's a / b / c / pri rule) / score (arm C).

## Decision rules (tier 1)
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

## Reported, not decided
- **Tier 2** = the stage's flagged candidates that are not in tier 1 (sharing < half of their members with every tier-1 flag):
  - its false-flag families on P4c;
  - the stage arm's own D right among the attainable reads, and its false moves on P4d.
- **Per deleted copy:** caught by tier 1, by tier 2 only, or by neither.
- **Also:** a per-family table, and the runtime.
