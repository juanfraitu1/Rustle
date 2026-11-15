# PREREG r1631: tie groups, genome-wide, on a sample the rules never saw (KB3781) + an assembly bar (2026-11-11)

Follows r1630 (ccm locus rule PASS on A119b and OR6737, both seen, both panel-scoped, grouping bars only). Three gaps are
tested here; written before any run. The recipe is frozen: `build PREFIX as98 0 1 RMSK:0.5 stranded ccm`, SEDEF co-duplication
in the scorer. Nothing is tuned on this sample. The layer sits in read assignment, which the thesis advisor has not approved:
opt-in only, not for `advisor_share`.

## Sample and input
KB3781 (`GCA_029281585.2_flnc_mm.bam`, default run `audit/g2/kb_r76`). No tie-group row (1603-1630) used it.
**Genome-wide input**: every molecule with `second_as >= 0.98 x best_as` (`best_as > 0`) in `kb.molecules.tsv`, all its
records (one `samtools view -N` pass). r1630 took only molecules over the panel; this is the product form.

## Part T: the grouping transfers
Truth as r1611 (gorilla RefSeq, protein-hit map `r1611/ens_map.tsv`), panel from `family_panel.py` on KB3781. No StringTie
or FLAIR run exists for KB3781, so M = **every** expressed panel member whose `hard_signal.py` signal is `indistinguishable`
(r1611 took only those both tools miss; this M is a superset in kind, fixed here).
- T1 grouped reach (B1g) >= 0.80 of M (r1603's 58/72).
- T2 purity counting co-duplicated loci (B2s) >= 0.90.
- T3 largest group (B5r): 1 chromosome, or >= 0.90 of its annotated loci in one family.
- T4 swallowed loci (D2) <= 0.05.
- T5 build <= 1 h; peak memory reported (no bar).

## Part A: is it an assembly improvement?
`bench/famcmp/tie_chains.py`: multi-exon tie-group models whose intron chain (contig, strand) is in neither `kb.gtf` nor
`kb.final.gtf` = **new chains**. Correct = exact RefSeq chain (`GGO.GCF_029281585.2_RefSeq.gtf.gz`). The layer only adds, so
no chain is lost.
- A1 new correct chains >= 15 (the r1627 yield bar on this sample).
- A2 share correct among new chains >= 0.070 (KB3781 default 2-read chains, r1627).
Reported, no bar: the same two numbers for models in groups of >= 2 loci; tie-group chains already in the default.

## Reading
- T and A pass: wire `tie-groups` as an opt-in stage (named, or `all --tie-groups`; needs an RMSK BED), then A119b genome-wide.
- T passes, A fails: the layer is a grouping product (which copies a read set cannot tell apart), not an assembly
  improvement. Say so; do not wire it as one; ask the user.
- T fails: the recipe does not transfer; stop.
- A119b cost: the panel build took 3,324 s. If the genome-wide A119b build passes 1 h, the stage is documented as
  gorilla-scale / panel-scale in Python and a Rust port is the next row, not a looser bar.

## Result (2026-11-11, job 4466351, `r1631/`): FAIL as registered
228,505 tied molecules, 1,928,348 records, 4,041 loci after repeats set aside (2,654), 2,623 groups; build 277 s, peak 1.9 GB.
- T1 9 of 9 PASS (M is only 9: 323 of 344 members are distinguishable on this sample). T2 148/148 = 1.000 PASS. T4 0.004 PASS. T5 PASS.
- T3 FAIL: the largest group has 43 loci on 9 chromosomes and none of its 27 annotated loci in one family.
- A1 15 new correct chains PASS (at the bar); A2 15 of 1,232 = 0.012 FAIL (bar 0.070). Groups of >= 2 loci: 7 of 499 = 0.014.
  Of the 1,864 tie-group chains, 632 are already in the default (333 correct).
Reading: the recipe's grouping does not transfer genome-wide (T3), and its models are not an assembly improvement (A2).
Not wired into the pipeline; no A119b genome-wide run.
