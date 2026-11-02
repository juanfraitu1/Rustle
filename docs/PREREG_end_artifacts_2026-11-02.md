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

## RESULT step 1 (2026-11-02): RI passes, chosen for the port
Baseline final.gtf: 211,392 models, = 39,277, chain 22.0 / 19.1.
| arm | dropped | = | chain sens / prec | multi/single |
|---|---|---|---|---|
| R (RTS, non-FSM) | 27,682 | 39,274 (-0.008%) | 22.0 / 22.1 | 1.01 |
| I (A-run, non-FSM, reads <= 2) | 14,903 | 39,277 | 22.0 / 20.6 | 1.03 |
| RI | 40,479 | 39,274 (-0.008%) | 22.0 / 23.8 | 1.00 |
All three pass; RI has the highest precision. Note the FSM exemption uses a reference: the port replaces it (step 2 design below).

## Amendment 1 (2026-11-02, before any step-2 code or numbers): the shipped rule is reference-free
Step 1's FSM exemption needs an annotation; the port replaces it by read support. Rule E (one rule, no tuning):
drop a final.gtf model if [some junction is RTS AND reads <= 2] OR [A share of the 20 bp downstream of the 3' end > 60% AND reads
<= 2]. RTS exactly as SQANTI3 6.1 rt_switching.py as qc calls it (min_match 8, allow_mismatch on (-a), wiggle 1, PATSEQLEN 10,
all junctions); A share as classification_classifiers.py (window 20, strand-aware).
Step-2 bars: (a) per-model detector agreement with job 4410027's RTS_stage and perc_A_downstream_TTS > 60 flags >= 95% each;
(b) rule E on A119b gw.LS final.gtf: chain precision >= +0.5, `=` loss <= 0.5%, multicopy/single removal <= 1.5x.
