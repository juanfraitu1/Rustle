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
