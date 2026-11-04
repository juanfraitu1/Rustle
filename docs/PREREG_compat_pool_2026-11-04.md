# PREREG r1526: compatible-chain admission for family members with reads but no two identical chains

User's goal text (2026-11-04), verbatim:
> Why: on the r1521 panel (A119b, strand vote shipped, base = fusedgw/R1), 107 of the 222 missed expressed family members have
> >= 2 local spliced reads but every read carries a different chain, so none reaches the pass-1 floor. The duplicate levers
> (ZMW, sequence) are closed. ChainAdmission::Compat exists but is not on the default path.
> Step 0 (before any run, decides whether to run at all): count, on the 107, how many have >= 2 reads whose chains are mutually
> compatible (one contains the other's introns, or they differ only in terminal exons). If fewer than 25, stop and record the
> ceiling; no arm.
> Arm: one genome-wide assemble with compatible-chain admission behind an env switch (RUSTLE_CHAIN_ADMIT=compat), reads admitted
> to a chain only when compatible AND on the same strand; the emitted chain is the longest member. Wall <= 45 min. Score with
> family_panel.py (definition T) and the chain-identity script.
> Bars (all must hold): complete families >= 209 + 10 (base 209); fused families <= 49; ordinary-member hit rate >= 0.829;
> base spliced chains kept >= 0.99; new chains' exact CAT match >= 0.05 (base 2-read chains 0.051); of the 107, >= 25 gain a
> same-strand model.
> PASS = ship as the assemble default; FAIL on any bar = stays opt-in, row closed.

## Step 0 (measured before this file, decides the arm)
Class = assemblable.tsv rows with ours = 0, maxP_chain < 2, maxT_chain < 2, n_local_spliced_P >= 1: **97 members** (the .out's
107 counted T-chain reads too; the 97 is the set scored). Same-strand spliced primaries over each gene, introns = N operations:
- goal's rule (>= 2 reads, one chain a contiguous sub-run of the other): **40 of 97** -> the arm runs (>= 25).
- the existing `--compat-support` rule (every junction of some read chain carried by >= 2 compatible reads): 23 of 97.
- 24 members have a single same-strand spliced read; 23 have >= 2 chains identical after 10-bp junction rounding
  (`--read-junction-snap`, not this row).
Scripts: scratchpad step0.py / step0b.py (copied to ~/r1526/).

## Fixed before any number
- The goal's mechanism (compatible reads pooled into the longest chain, admitted at the floor, same strand) IS the existing
  `--direct-transfrags-union` (register 1266, 2026-10-11: every default chain plus every maximal transfrag with the reads of the
  transfrags it contains; not adopted then on gffcompare precision, never scored on the family panel). No new code; the env
  switch is the existing `RUSTLE_DIRECT_TRANSFRAGS=union`. r1266's -0.8 to -1.5 chain precision is the prior for the CAT-match bar.
- Base = fusedgw/R1/x.gtf (strand vote on, layers off; chains identical to A). Arm = fusedgw/U/x.gtf: one genome-wide
  `rustle_pipeline.sh assemble`, RUSTLE_ASSEMBLY_LAYERS=off RUSTLE_DIRECT_TRANSFRAGS=union, strand vote default, current main
  binary (~/rustle_target_test/release, d32c6900). Wall bar 45 min on the assemble step.
- Scores: `bench/famcmp/family_panel.py` definition T (complete, fused, ordinary = hit_rate_mapq0_share<0.5);
  `bench/famcmp/chain_cat.py BASE ARM CAT` (base spliced chains kept; new chains' exact CAT Liftoff chain match, the r1524 proxy);
  `bench/famcmp/score34.py --one-chain` (of the 97, same-strand model overlapping the gene's exons).
- Bars as in the goal, with "of the 107" read as of the 97.
