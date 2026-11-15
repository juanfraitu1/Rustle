# PREREG r1632: sibling-projected models as the guide of `assemble --guide` (KB3781, 2026-11-11)

User's choice after r1631 ("Sibling models as guides"). Written before the run.
Guide = `r1627/kb/projected.gtf` (every family model projected onto the other members, r1519 parameters, all tiers).
Run = `utilities guide --gtf kb.gtf --guide projected.gtf --bam KB3781` (rule G1d as shipped; nothing changed).
Prior from the record, stated before running: the chains this can add are r1627's classes S + T (a projected chain with a
full-chain primary record), 16 correct of 263 = 0.061, so N2 is expected to fail unless the guide's removals change it.

Scoring (`bench/famcmp/guided_sibling.py`), RefSeq-exact chains, only chains whose span overlaps a projected model on its
contig and strand (the guide covers family loci only; the mode's "outside the guide, < 3 reads" drop elsewhere is an
artefact of a partial guide and is reported, not scored):
- N1 new correct chains (guided, not default) >= 15.
- N2 share correct among new chains >= 0.070 (KB3781 2-read default chains).
- L1 new correct minus lost correct > 0, and share correct among lost chains (default, not guided) <= 0.070.
All three pass = PASS. Otherwise FAIL: the family-sharing direction is closed for assembly; ask the user.

## Result (2026-11-11, job 4466383, `r1632/`): FAIL as registered
Guide verdicts on 72,460 de novo models: Keep 56,370, NovelLow 14,574, NewIntronLow 1,283, Fragment 214, Retained 19; 263 guide transcripts added (175 s).
- N1 16 new correct PASS. N2 16 of 263 = 0.061 FAIL. L1 FAIL: 1,631 chains lost at family loci, 162 of them correct (0.099), net -146.
