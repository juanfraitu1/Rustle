# PREREG row 1724: full-length exemption from the isoform-fraction rule (2026-11-14)

User go 2026-11-14 ("yes lets test the fraction rule exemption"). Written before any figure of this row was read.

## Question
The fraction rule (`--polish-isoform-fraction 0.02`) accounts for about 40-45% of the lost annotated minor chains
(`docs/PREREG_alt_rules_2026-11-14.md`). Are the chains it drops real when their reads are full-length at both ends?

## Design
- Samples: A119b (discovery) and HG002 on CHM13 (held out, decides). Human libraries only: the cap gate refuses both gorilla
  libraries (row 1720).
- One job assembles each sample twice with the current binary: the shipped driver, and a copy of today's `tools/` whose
  only change is `--polish-isoform-fraction 0.0`. Checks before scoring: the fraction-off `PREFIX.assemble.polish` shows 0.0;
  the rerun default chain set equals the stored one (`r1617/default/a.gtf`, `r1692_trunc5_hg/hg.gtf`), else earlier
  baselines are not comparable and that is reported.
- Candidates: multi-exon chains in the fraction-off `PREFIX.gtf` that are absent from the default `PREFIX.gtf`. Default
  chains that go missing with the rule off are reported. "Default plus exempted chains" approximates a real exemption:
  inside the binary other polish rules would see a different chain set.
- Signature (frozen, existing code, `bench/famcmp/capg.py exempt`): 5' = a capped cluster of the chain's own reads
  (`clusters`: >= 3 capped reads within 10 bp, >= half of the informative reads there); 3' = the `earlyends` rule (commonest
  read end holds >= half of the reads within 10 bp, polyA signal upstream, no A-rich run after). Chains with no informative
  5' read (G upstream) are their own group, "undecidable".
- Decision set E: candidates with both ends that are NOT a terminal sub-run of a default chain (those are what `capstarts`
  and `earlyends` already reach). Truth: RefSeq chain identity.

## Bars (HG002; all four for PASS)
- B1 quality: annotated share of E >= the sample's default chain precision.
- B2 depth control: annotated share of E >= 1.5 x that of the same number of top-read candidates without the signature
  (decidable, not a sub-run). Without it a pass could only mean "more reads".
- B3 yield: annotated chains in E >= 1% of the default's annotated chains.
- B4 combined precision: (default + E) >= default - 0.005.

Secondary, reported, not deciding: the 3'-only signature (any 5' state, not a sub-run) on all four samples, gorilla from
the stored `r1614/{default,nofrac}/or.gtf`; A119b figures.

Measurement only: no driver stage, flag or default change follows from this row without a user ruling.

## Outcome
OUTCOME (job 4471916) FAIL, A119b / HG002: candidates 58,118 at 0.051 / 81,483 at 0.056 annotated; decision set 4,625 at 0.075 / 10,016 at 0.071 (bars 0.176 / 0.234); equal-count top-read no-signature control 0.089 / 0.090, so the signature does not beat read depth; annotated chains 345 (bar 387) / 710 (bar 357); combined precision 0.174 / 0.224 (HG002 bar 0.229). Rerun default equals the stored one on both. Secondary 3'-only arm 0.055 / 0.059 / OR6737 0.059. The fraction rule drops chains that are about 5% annotated whatever their ends show; no exemption.

All groups sit between 0.02 and 0.09 annotated on both samples (both ends 0.075 / 0.071, cap only 0.069 / 0.075, polyA end only 0.053 / 0.056, neither 0.044 / 0.046). Only 7 candidates per sample are terminal sub-runs of a default chain, so this set does not overlap what `capstarts` / `earlyends` reach. With the rule off 290 / 584 default chains go missing (27 / 43 annotated). Tables: workspace `r1724_fracexempt/{a,hg,or}.txt`. Truth is RefSeq chain identity, which understates unannotated real isoforms; the control uses the same truth, so B2 does not depend on that.

# PREREG row 1731: are the ends of the exempted chains real? (CAGE peaks, polyA sites)

User go 2026-11-14 ("ok lets preregister and run"). Written after row 1724's outcome and before any figure of this row.

Row 1724 scored chains by RefSeq chain identity, which understates real unannotated isoforms. This row scores the ENDS of
the same candidates by two independent yardsticks already on disk (both human, on CHM13): refTSS CAGE peaks
(`r1716_sqanti/refTSS.hs1.bed`, 5' position inside a peak; background 100 bp downstream) and GENCODE V44 polyA sites
(`r1719_polya/sites.hs1.bed`, within 25 bp of the 3' end; background 200 bp upstream).

- Same candidates and signature as row 1724 (`capg.py exempt` rerun with two added columns: 5' position = the capped
  cluster, else the commonest read start; 3' position = the commonest read end), then `capg.py exempt-ends`.
- Groups, terminal sub-runs excluded: E unannotated; E annotated (the reference rate for real chains with this signature);
  control = the same number of top-read unannotated candidates without the signature.
- Bars on HG002 (A119b reported): B1 5' share of unannotated E >= half that of annotated E; B2 the same for the 3' end;
  B3 share with both ends supported >= 2 x the control's.
- Reading fixed in advance: a PASS says the ends are real, not that the intron chains are. No exemption, stage or default
  follows from this row; it only qualifies row 1724's 7% figure.

## Outcome (row 1731)
OUTCOME (job 4472005) FAIL as registered on HG002 (B3), A119b / HG002: unannotated decision set 4,280 / 9,306; 5' in a CAGE peak 0.648 / 0.907 (annotated members 0.788 / 0.879; background 0.07 / 0.05); 3' at a polyA site 0.749 / 0.838 (annotated 0.843 / 0.854; background 0.01); both ends 0.497 / 0.762 vs control 0.237 / 0.408 (ratio 2.10 PASS / 1.87, bar 2). The ends of these chains are as well supported as those of annotated chains, on HG002 nearly identical; the HG002 control is also high (5' 0.799). Row 1724 tables reproduced exactly; its control row moves by 4 chains on rerun (ties at the read-count cut, row order). Ends real, chains unproven: no exemption. Tables: workspace `r1724_fracexempt/{a,hg}.ends.txt`.

# PREREG row 1733: are the intron chains of the exempted chains real? (known junctions, an independent sample)

User go 2026-11-14 ("lets test for that and add the new row"). Written after rows 1724 and 1731, before any figure of
this row. No short-read junction set is on disk, so two chain-level yardsticks from what we have:

- Known junctions: every intron of the chain is an intron of some RefSeq transcript (a new combination of known junctions).
- Replication: the exact intron chain is spelled by >= 1 primary read of the OTHER human sample (A119b testis for HG002
  candidates and the reverse; different individual, tissue and library, both on CHM13).

Same candidates and signature as row 1724 (`capg.py exempt ... OTHER_BAM`, then `capg.py exempt-junc`); terminal sub-runs
excluded. Groups: E unannotated; E annotated (reference: what real chains with this signature show; replication only);
control = unannotated no-signature candidates matched to E unannotated by read count (seed 0), because replication grows
with depth and row 1724's top-read control would be biased upward.

Bars on HG002 (A119b reported; all three for PASS): B1 known-junction share of unannotated E >= 1.5 x the control's;
B2 replicated share >= half that of annotated E; B3 replicated share >= 1.5 x the control's.

Limits stated in advance: tissue differs, so absence in the other sample is weak evidence against a chain; alignment
artifacts at repeats can replicate. A PASS makes the chains credible; whether to exempt them at row 1724's precision cost
(RefSeq terms) is then the user's call. No stage or default follows from this row.

## Outcome (row 1733)
OUTCOME (job 4472019) FAIL, A119b / HG002: unannotated decision set 4,280 / 9,306; all introns annotated 0.074 / 0.068 vs read-count-matched control 0.142 / 0.095 (below it); chain spelled by a read of the other sample 0.318 / 0.374 vs annotated members 0.710 / 0.599 and control 0.341 / 0.394 (not above it). About 93% of these chains carry an intron RefSeq does not have, and they replicate no better than fraction-rule drops without the signature. With row 1731: real ends, intron chains no better than the rest of what the rule drops. No exemption; closed. The control has 4,174 / 8,985 chains (some read-count strata ran short). Tables: workspace `r1724_fracexempt/{a,hg}.junc.txt`.
