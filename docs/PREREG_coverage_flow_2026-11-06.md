# PREREG: coverage-conserving flow for tie resolution (row 1554, 2026-11-06 22:15, user-requested)

## Question
Read-star assigns a read when one copy wins by a margin; ties and ambiguous reads stay unassigned. A bipartite
read-copy graph with per-copy capacity from the copies' own depth (min-cost flow; equivalently EM with copy
weights) can resolve ties toward copies the unique evidence supports. Does it agree with read-star where
read-star is sure, and how many ties does it resolve? Answer lever only: no alignment is saved.

## Set-up (MCL1; row 1544 ctl_genomic.paf = minimap2 -x splice, 15,828 reads x 31 genomic windows, -N 31)
Per read and copy: AS = best alignment score of that read on that copy's window (missing = no alignment).
EM: w_c uniform; E-step P(c | r) proportional to w_c * exp((AS_rc - AS_r,best) / T) with T = 10 AS units (the
scale of read-star's own margin rule, as_margin); M-step w_c = mean_r P(c | r); 50 iterations. A read is
"flow-assigned" to c iff P(c | r) >= 0.9.
## Bars
- Agreement: of the baseline-assigned reads present in the PAF, >= 90% flow-assigned to the same copy and
  <= 2% to a different copy.
- Ties: report how many of the baseline-tied reads become flow-assigned, and the distribution of their copies
  (concentration on a few copies = the rich-get-richer artifact of coverage priors). No truth is available on
  A119b ties; a PASS here is "consistent with read-star", not "correct".
## Prediction
Agreement passes (large margins dominate the prior). Ties: 30-70% resolved, concentrated on the 3-5 copies with
the most unique reads. The row cannot establish correctness; if agreement passes, the next step would be the
same procedure on KB3781, where haplotype truth exists (memory: O2 ms-tie scope), before any adoption.
