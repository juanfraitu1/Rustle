# PREREG END ARTIFACTS (register 1500, 2026-11-02, before computing any arm)

Why: r1499 (SQANTI3 default rules) raised chain precision 19.1 -> 25.4 but lost 6.4% of `=`. SQANTI3 is used here only to choose a
rule; the shipped rule will be Rust code of ours (no SQANTI3 dependency).
Step 1 (oracle, existing flags of job 4410027 on A119b gw.LS final.gtf): drop from final.gtf
- R: models with RTS_stage TRUE that are not full-splice_match;
- I: models with perc_A_downstream_TTS > 60 that are not full-splice_match AND have reads <= 2 (final.gtf `reads`);
- RI: R or I.
Bars per arm, gffcompare vs chains/ref.gtf: (1) chain precision >= +0.5 points, (2) `=` loss <= 0.5%, (3) multicopy (SD90 overlap)
removal share <= 1.5x single-copy. The passing arm with the highest chain precision is ported (step 2); none pass -> FAIL, no port.
Step 2 (own code): RTS = a junction whose donor/acceptor flanks hold a direct repeat (SQANTI's definition: >= 8 bp match within
the 8 bp windows... ported exactly from its source) and intrapriming = A share of the 20 bp after the 3' end, as a final.gtf layer;
pass = per-model agreement with the SQANTI3 flags >= 95% and the same three bars on A119b genome-wide.
