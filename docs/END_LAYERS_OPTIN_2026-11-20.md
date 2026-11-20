# Row 1846: the three end layers of `final.gtf` become opt-in (user decision, 2026-11-20)

User (2026-11-20): "lets turn all those defaults you mentioned into opt ins" (the single-exon model layer, the 3' UTR
extension and the median 5' end). Basis: row 1831 (`docs/PREREG_mono_end_rules_2026-11-20.md`): mono-model pays on the two ape
libraries and fails on zebrafish, mouse ES and heart; the 3' extension has an effect on ape testis only; the median 5' end
helps A119b and mouse ES and hurts the other libraries; row 1834: no cap-signal gate predicts where it helps.
Change: `tools/rustle_pipeline.sh`, new `RUSTLE_FINAL_END_LAYERS` (`tssend`, `tesext`, `monomodel`, `all`; unset = none).
`RUSTLE_FINAL_END_LAYERS=all` = the output before. The sub-chain layer of `final.gtf` (`RUSTLE_FINAL_SUBCHAIN`), the chain
score filter and the 5' truncation rule are untouched. `PREFIX.gtf` is not affected by construction. Workspace
`r1846_end_layers_optin`.

## Check written before the runs (not a bar for the decision, which is made)
Row 1831 switched one layer off at a time, on seven datasets, before the junction replacement was a default. Runs: the new
default and `RUSTLE_FINAL_END_LAYERS=all` on the twelve datasets of rows 1837-1845. Reported for `x.final.gtf` with row 1831's
scorer: multi-exon chains and correct chains, single-exon models / precision / truth found, 5' and 3' ends within 50 bp among
correct chains, exonic base sensitivity and precision. Parity: the `all` arm equals row 1845's `jr20` byte for byte; `x.gtf`
is identical between the arms.
Flag to the user: a dataset where the new default loses more than 0.5% of the correct chains or more than 1 point of chain
precision in `final.gtf`; and the size of what the ape libraries give up (ends, bases, single-exon), which row 1831 predicts.
Known interaction: the chain score filter and the truncation rule read end features; row 1831 saw 3-49 MORE correct chains
with the median 5' end off.
