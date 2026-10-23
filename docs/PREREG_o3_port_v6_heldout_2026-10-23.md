# Pre-registration: the O3 port v6 on a fresh gorilla fibroblast held-out (P5), 2026-10-23

**Status:** committed before the 20 families are drawn; no data of these families has been aligned or scored.

## Why
- **Development.** Port v6 passes both development panels (`docs/PREREG_o3_port_2026-10-23.md`, development log):
  - P1: D right 6,129 = 83.9% of C; false moves 0.02%.
  - P2: 0 of 53 families with a false flag.
- **Held-out so far.** P3 was spent on v4, and every change since was made on P1 / P2, so v6 has no held-out evidence.
- **Correlation clustering** passes development too (P1 5,949, P2 1 / 53) and is better in the one-copy simulation when the sample is
  another individual.
- **User constraints** (2026-10-23): "lets stay in gorilla, try the held out". Earlier: keep it small ("this might be overly
  complicating").

## Panel (fixed before drawing)
- **Families:** 20 drawn with `random.Random(1).sample` from P4's registered draw (`/mnt/linuxdisk/tmp/rna_allele/p4/panel.json`).
  P4's draw has 60 two-copy fibroblast families in neither the 53 development families nor refabsent's 34, and was never aligned
  or scored (`docs/PREREG_o3_two_tier_p4_2026-10-23.md`, "Not run").
- **Deleted copy:** P4's (Amendment 7's rule).
- **Reads:** P4's scored reads of the 20 (`link_test.py reads`: <= 500 primaries per copy, seed 1, the fibroblast BAM).
- **Deletion arm:** the `_pri` with only these 20 copies masked.
- **Control arm:** the unmasked `_pri` with every copy.
- **Alignment command:** Amendment 7 / 9's.

## Code (frozen)
- **Port:** `bench/o3_graph/o3port.py` as of e11cf8c6.
  - **Primary arm:** v6, the default options (label propagation).
  - **Secondary arm:** `--clusterer cc`.
  - Both use the stage's nets and products of the same arm (IsoCon side).
- **Stage:** `candidate_copies` frozen at `p4/bin/candidate_copies` (sha1 2e3cd592). It reproduces Amendment 15's products byte for
  byte on development batch g0. Batches of 10 families.
- **Scorers:**
  - Deletion arm: `merge_test.py score`, arm M.
  - Control arm: `control_test.py` lift / classify, Amendment 9's a / b / c / pri rule.

## Decision rules (the primary arm)
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

## Reported, not decided
- Look-harder flags by class (control) and by source (deletion).
- The landing loci.
- Each deleted copy's outcome: certified, look harder, in a slot, or no pile.
- Runtime.
