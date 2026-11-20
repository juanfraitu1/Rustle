# Row 1845: junction replacement at 20 x support becomes the driver default (user decision, 2026-11-20)

User (2026-11-20): "lets make the 20x junction replacement the default and update the results". Basis: row 1819
(`docs/PREREG_junction_replace_2026-11-19.md`: at ratio 20 at least 0.995 of the correct chains kept and higher precision on
seven check datasets, measured on `PREFIX.gtf`) and row 1835 (ported as `--polish-junction-ratio`, parity with row 1819).
Change: `tools/rustle_pipeline.sh`, `RUSTLE_POLISH_JUNCTION_RATIO` unset = 20; `RUSTLE_POLISH_JUNCTION_RATIO=0` = the command
before. `copy_assign`'s own default stays 0. Workspace `r1845_junction_ratio_default`.

## Check written before the runs (not a bar for the decision, which is made)
Rows 1819 / 1835 read `PREFIX.gtf`; `final.gtf` (the chain score filter, fitted without the replacement) was not measured
with it. Runs: the new default on the twelve datasets of rows 1837-1843, against their `base` (row 1843, same binary source).
Reported per dataset for `x.gtf` and `x.final.gtf`: correct chains, precision, correct lost and wrong removed; SIRV contigs of
the mouse ES library apart (`r1801_sirv_wrong/sirv.truth.gtf`). Parity: `RUSTLE_POLISH_JUNCTION_RATIO=0` reproduces the old
`base` byte for byte (fly, mouse ES).
Flag to the user if, in `final.gtf` on any dataset, more than 0.5% of the correct chains are lost or precision falls.
Expectation: `x.gtf` as row 1819; in `final.gtf` the same direction, smaller, since the filter already removes part of the
shifted-site chains.
