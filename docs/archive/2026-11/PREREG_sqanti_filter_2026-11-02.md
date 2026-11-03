# PREREG SQANTI3 RULES FILTER (register 1499, 2026-11-02, before the run)

Why: r1498 cut the contained-chain share of `PREFIX.final.gtf`; SQANTI3 flags other artifact classes (intrapriming, RT-switching,
non-canonical junctions) from the genome and reads, so a rules filter could be a fifth final.gtf layer for any species.
Input: A119b genome-wide `finalsub/gw.LS/x.final.gtf` (today's default final.gtf, 211,392 transcripts); reference for SQANTI3
classification and gffcompare = `chains/ref.gtf` (CHM13 RefSeq, as r1498); genome chm13v2.0.fa. chr20 smoke first (wall check).
Arms: F0 = final.gtf as is; F1 = `sqanti3_filter.py rules` (default rules, as r1482) output on the SQANTI3-corrected models of F0.
Report (no decision): SQANTI3 categories and artifact flags of F0, split base (model also in PREFIX.gtf, same id and exons) vs layer,
and single-copy vs `multicopy "true"`.
Decision: F1 becomes a proposed opt-in final.gtf layer if, vs RefSeq by gffcompare, (1) intron-chain precision rises >= 0.5 points,
(2) `=` falls <= 0.5%, AND (3) the share of multicopy models removed is <= 1.5x the share of single-copy models removed (the filter
must not target family copies). Else FAIL. The filter's own artifact share is not a bar (it removes what it flags).

## Amendment 1 (2026-11-02, after the chr20 smoke, before any genome-wide output)
The smoke showed `multicopy` is "false" on every assemble-stage model (copies are adjudicated later), so bar (3) was unmeasurable.
Multicopy = any exon overlapping `CHM13_SDs.SD90.bed.gz` (segmental duplications >= 90%). Base = same exon chain in PREFIX.gtf,
any id (the layer run renames ids). gffcompare tmap path fixed (written beside the input). Bars and arms unchanged.

## RESULT r1499 (2026-11-02, job 4410027): ❌ FAIL (bar 2)
| arm | transcripts | = | c | chain sens / prec |
|---|---|---|---|---|
| F0 final.gtf | 211,392 | 39,277 | 11,421 | 22.0 / 19.1 |
| F1 rules filter | 151,317 | 36,779 (-6.4%) | 8,974 | 20.6 / 25.4 |
Bar 1 pass (+6.3 points), bar 2 FAIL (= -6.4% vs <= 0.5%), bar 3 pass (multi/single removal 1.09).
Removals: RT-switching 26k, intrapriming (perc_A_downstream_TTS > 59) 37k; every non-FSM removal also cites "NA value in min_cov"
(no short-read junction coverage given). Base models flagged most (45.8% removed vs 20.8% of layer-changed models; intrapriming
26.8% vs 10.8%): the 3' end layer already moves many A-run ends. Not adopted. A follow-up would need its own prereg (e.g. RTS only,
or intrapriming with an FSM/read-support exemption).
