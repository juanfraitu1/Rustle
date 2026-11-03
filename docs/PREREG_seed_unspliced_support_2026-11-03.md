# PREREG UNSPLICED SECONDARIES AS SUPPORT ONLY (register 1520, 2026-11-03)

**Written before any number of the arm.** r1518 (`docs/PREREG_seed_unspliced_primary_2026-11-03.md`) dropped unspliced secondaries
from the seeding pool: every array-wide single-exon blob went away and families were untouched, but CAT single-exon genes hit fell
757 -> 342 (Amendment 1 rule), because the same near-tied placements lift single-exon pseudogene copies over the pass-1 floor.
User (2026-11-03): "can we register and run the next step then?"

## Change (`RUSTLE_SEED_UNSPLICED=support`, opt-in, byte-identical unset; `primary` = r1518's rule)
The streaming reader holds unspliced SECONDARY records apart (`Pass1Acc::unspliced_sec`). Pass 1 clusters the PRIMARY unspliced reads
by span single-linkage as before (`cluster_unspliced`, floor 1), then each held secondary adds one read of support to the first
primary cluster whose span it overlaps and is otherwise discarded; clusters below `min_reads` after that are dropped
(`cluster_unspliced_supported`). Boundaries come from primaries only. A secondary can therefore lift a copy over the floor
(r1060's fill-in) but can neither seed a locus with no primary nor chain two primary clusters into one blob. Spliced
secondaries unchanged. Unit test: 1 primary + 3 overlapping secondaries -> one skeleton n=4 with the primary's span; 3 secondaries
with no primary -> nothing.

## Arms and bars: r1518's, unchanged (A = fusedgw/A; C = this rule; panel chr4+chr20 assemble first, then genome-wide)
1. Single-exon models >= 20 kb with a zero-primary run >= 3 kb: 0 (A 61; cannot be otherwise if the mechanism is as stated).
2. chr20 SQANTI single-exon: artifact classes <= 253; FSM >= 0.98 x 16 (i.e. >= 16).
3. Pair precision >= 0.7158, recall >= 0.7307, fused loci <= 547.
4. CAT single-exon genes hit (Amendment 1 rule, model span <= 2x gene) >= 0.98 x 757 = 742. Reported beside: the r1518 B value 342.
5. Spliced models identical to A >= 0.99.
PASS = all five -> default (`--seed-pool good` gets the support rule). Any failed = FAIL, opt-in.

## RESULT (2026-11-03 21:47): FAIL as registered (bar 2 FSM half, bar 4); the mechanism holds and the cost is a third of r1518's
| bar | A | C (support) | B (r1518 drop) | verdict |
|---|---|---|---|---|
| 1 zero-primary-run models | 61 | 0 | 0 | PASS (USP17L span: 2 local models, was 1 blob / 0) |
| 2 chr20 single-exon artifacts / FSM | 253 / 16 | 233 / 15 | 227 / 14 | artifact PASS, FSM FAIL by one |
| 3 pairs precision / recall / fused | guards 0.7158 / 0.7307 / 547 | 0.726 / 0.735 / 545 | 0.728 / 0.742 / 545 | PASS |
| 4 CAT single-exon genes hit (span <= 2x) | 757 | 658 (0.869) | 342 (0.452) | FAIL |
| 5 spliced models identical | 219,503 | 219,503 | 219,503 | PASS |
Single-exon models 16,501 -> 15,316 (B 14,498). Reading: letting the tied secondaries support but not seed keeps two thirds of the
pseudogene fill-in r1518 lost (658 of the 757 - 342 = 415 gap recovered, 76%) while every blob stays gone. The remaining 99 lost
genes are copies with < 3 primaries whose cluster the secondaries used to CREATE, i.e. the r1060 class itself. Stays OPT-IN
(`RUSTLE_SEED_UNSPLICED=support`); a default change needs a bar on bar 4 that the user accepts as the price of the blobs.
