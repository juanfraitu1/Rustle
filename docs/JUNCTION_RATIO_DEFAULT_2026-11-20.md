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

## Result (2026-11-20): the default holds in `final.gtf` on all twelve datasets; nothing to flag

Jobs 4490621-4490632, 24 runs exit 0, shared binaries `~/rustle_target` rebuilt from `3e26f1e7`;
`r1845_junction_ratio_default/cmp.NAME.txt` (`cmp.py`). Parity: `RUSTLE_POLISH_JUNCTION_RATIO=0` gives `x.gtf` and
`x.final.gtf` byte-identical to row 1843's `base` on all twelve.

`x.final.gtf`, ratio 0 -> ratio 20 (the new default): correct chains / precision, correct lost, wrong chains removed.
| dataset | before | after | correct lost | wrong removed |
|---|---|---|---|---|
| fly simulated | 7,041 / 0.9643 | 7,016 / 0.9669 | 0.36% | 8.0% |
| mouse simulated | 27,159 / 0.9534 | 27,041 / 0.9587 | 0.43% | 12.3% |
| zebrafish simulated | 12,381 / 0.8910 | 12,360 / 0.8949 | 0.17% | 4.1% |
| mouse ES | 17,730 / 0.3472 | 17,678 / 0.3669 | 0.29% | 8.5% |
| mouse heart | 13,091 / 0.5768 | 13,068 / 0.5899 | 0.18% | 5.4% |
| mouse adrenal | 14,394 / 0.5779 | 14,361 / 0.5912 | 0.23% | 5.5% |
| mouse cortex | 15,679 / 0.4868 | 15,655 / 0.5002 | 0.15% | 5.4% |
| OR6737 | 25,635 / 0.3775 | 25,565 / 0.3835 | 0.27% | 2.8% |
| A119b | 38,091 / 0.2187 | 37,936 / 0.2207 | 0.41% | 1.6% |
| chimpanzee | 20,847 / 0.4155 | 20,807 / 0.4195 | 0.19% | 1.8% |
| bonobo | 17,884 / 0.3367 | 17,840 / 0.3400 | 0.25% | 1.7% |
| KB3781 | 25,745 / 0.3733 | 25,667 / 0.3802 | 0.30% | 3.2% |
`x.gtf` moves alike (correct lost 0.15-0.42%, wrong removed 1.2-7.8%). No dataset loses more than 0.5% of its correct chains
and precision rises on every one. Expectation "smaller in `final.gtf`" was wrong: the share of wrong chains removed is as
large or larger after the filter.

**SIRV contigs of the mouse ES library** (all 22 contigs, classes of `r1801_sirv_wrong/klass.py`):
| | reported | correct of 60 | wrong | shifted site <= 10 bp | other new intron | piece | combination |
|---|---|---|---|---|---|---|---|
| `x.gtf`, ratio 0 | 90 | 53 | 37 | 25 | 5 | 6 | 1 |
| `x.gtf`, ratio 20 | 80 | 54 | 26 | 14 | 5 | 6 | 1 |
| `final.gtf`, ratio 0 | 76 | 48 | 28 | 22 | 5 | 0 | 1 |
| `final.gtf`, ratio 20 | 67 | 48 | 19 | 13 | 5 | 0 | 1 |
The default removes 9 of the 22 shifted-site false transcripts of `final.gtf` and no correct one (precision 0.632 -> 0.716).
Left: 13 shifted sites with more than a twentieth of the true site's support, and 5 with an intron in an unspliced spike-in
or elsewhere. A read never joins a SIRV contig to a mouse contig in the assembly (66,661 primary reads on the SIRV contigs, 1
with a supplementary part on a mouse contig; supplementary records are not assembled).
