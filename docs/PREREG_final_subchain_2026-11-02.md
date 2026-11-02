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
