# PREREG FINAL SUB-CHAIN LAYER (register 1498, 2026-11-02)

Why: genome-wide A119b vs native CHM13 RefSeq (gffcompare): ours `=` 38,754 vs StringTie 32,006, but contained chains `c` 16,251
(6.9% of 234,081) vs 3,131 (1.3%). `--polish-subchain drop` (r1130/r1138, opt-in) halved `c` on A119b held-out at <= 0.3% chains.
Change: `RUSTLE_FINAL_SUBCHAIN=drop` adds `--polish-subchain drop` to the layered run that writes PREFIX.final.gtf only
(PREFIX.gtf / families input unchanged by construction).
Run: A119b genome-wide `rustle_pipeline.sh assemble` (default good pool), arm L (layers, today's default) vs arm LS (layers + drop).
Metrics (gffcompare -r A119b.chm13v2.0_RefSeq.gtf, PREFIX.final.gtf): `c` share of query transcripts; `=` count; intron-chain
sensitivity/precision.
Decision: LS becomes the default (RUSTLE_FINAL_SUBCHAIN unset = drop) if (1) the `c` share falls >= 30% relative, (2) `=` falls by
<= 0.5%, and (3) intron-chain precision does not fall. Else FAIL, stays opt-in. Also checked: PREFIX.gtf byte-identical between arms.

## RESULT r1498 (2026-11-02): ✅ PASS (job 4407880)
| A119b PREFIX.final.gtf vs CHM13 RefSeq | query | `=` | `c` | `c` share | chain sens / prec |
|---|---|---|---|---|---|
| L (three layers) | 228,717 | 39,452 | 23,060 | 10.1% | 22.1 / 17.6 |
| LS (+ sub-chain drop) | 211,392 | 39,277 | 11,421 | 5.4% | 22.0 / 19.1 |
(1) `c` share -46% (bar -30%) ✓; (2) `=` -175 = -0.44% (bar <= 0.5%) ✓; (3) chain precision +1.5 ✓. PREFIX.gtf byte-identical
between arms (chr20 smoke and genome-wide). chr20 smoke: `c` 518 -> 274, `=` 1,097 -> 1,092.
Default: `RUSTLE_FINAL_SUBCHAIN` unset = drop (off restores the three-layer final.gtf). Note the three layers alone raise `c` (PREFIX.gtf
6.9% -> final 10.1%); with the drop final.gtf is 5.4% (StringTie 1.3%).
