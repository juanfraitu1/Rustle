# PREREG: read-overlap graph on the ambiguous pool (row 1553, 2026-11-06 21:50, user-requested)

## Question
Reads from one copy overlap each other with HiFi-level mismatch (~0.1%); reads from different copies at ~99%
identity differ ~1 per 100 bp. Connected components of a strict overlap graph should then be copies, and one
alignment per component would replace one per read. Does the graph separate MCL1's copies, or fuse?

## Set-up (MCL1, row 1544 reads: 15,828 of 30,034 molecules; baseline gw7: 35 assigned, 41 tied, 15,752 ambiguous)
`minimap2 -X -c --eqx -x map-hifi -t 16 reads.fa reads.fa`; edge iff overlap >= 500 bp and NM / aligned <= 0.002;
union-find components. Component -> copy by the baseline-assigned reads it contains.
## Bars
- Time: all-vs-all <= 160 s at 16 threads (the 469 s control / 3).
- Purity: no component holds baseline-assigned reads of two different copies (<= 1 violation).
- No fusion: largest component <= 20% of aligned reads (a giant component means copies fused).
- Coverage: the 35 assigned reads fall in components that together hold >= 50% of the reads (otherwise the
  components are too small to replace per-read alignment and the cost is not saved).
## Prediction
Time passes (ava on 16k HiFi reads is ~1-3 min). Purity is plausible at this strictness. Fusion vs fragmentation
is the coin flip, and the pre-registered expectation is FAIL: either a giant component (reads bridging copies
through regions without PSVs) or thousands of tiny components (strict rule + 1.7 kb reads with few overlaps per
copy). Row 1552 showed graph bubbles do not separate these copies; overlaps may do better because they use
base-level identity, which is why the arm is run.

## Result (2026-11-06 22:05): FAIL by fusion, as predicted; row closed
r1553 (job 4450724): all-vs-all 204 s at 16 threads (bar <= 160 FAIL), 2.93M overlaps, 498k edges under the strict
rule. Components 2,267 of which 2,225 singletons; the largest holds 11,524 of 15,828 reads (72.8%; bar <= 20%
FAIL). The 35 assigned reads sit in 10 components with 1 purity violation (bar <= 1, PASS but inside the giant
component). Coverage 73.2% PASS only through the giant component. Reading: reads bridge copies through stretches
without a PSV, so even a 0.2% mismatch rule percolates into one component; both failure modes appeared at once
(a giant component plus 2,225 singletons). Lever closed. No further rule tried (post hoc).
