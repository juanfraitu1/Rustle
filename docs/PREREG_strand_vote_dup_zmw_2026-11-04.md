# PREREG r1525: ship the strand vote; narrow the duplicate rule to different ZMWs

User's goal text (2026-11-04), verbatim:
> Mechanism: (1) single-exon models take their strand from the read-orientation vote (RUSTLE_READ_STRAND on, pre-fix binary
> da6c15ee behaviour); (2) coordinate-duplicate admission (RUSTLE_DUP_ADMIT) restricted to duplicates from different ZMWs,
> measured alone against base.
> Panel: r1521, 420 primates-level families, A119b, judged under T; the 34 assemblable misses (~/famcmp/assemblable.tsv);
> r1518 single-exon bars (bench/seed_unspliced_score.py on A vs arm, full output kept).
> Bars for (1): spliced chains identical >= 0.99; complete families >= 204; ordinary-member hit rate >= 0.819; CAT single-exon
> gene hits >= base (r1518 bar 4 as registered); the single-exon model count drop explained by class.
> Bars for (2): new chains' exact CAT chain match >= 0.05 (base's own 2-read chains); base chains kept >= 0.99; >= 15 of the 34
> gain a model.
> Wall: one genome-wide run per part.
> PASS (1) -> RUSTLE_READ_STRAND default on. PASS (2) -> default on. A part that fails stays opt-in; the fused-family cut is
> closed as de novo.

## Fixed before any number
- Base = fusedgw/A/x.gtf. Arm (1) = fusedgw/R1/x.gtf (job 4435820, r1524: pre-fix binary + RUSTLE_READ_STRAND=1, layers off);
  its chain identity (1.0000), panel (complete 209, ordinary 0.829) were seen in r1524 and are restated, not re-measured.
  New for (1): the full seed_unspliced_score.py output (bar 4 as registered: CAT single-exon genes hit, A vs R1) and the
  class breakdown of the 2,993 fewer single-exon models.
- Arm (2) = one genome-wide run, fixed binary with RUSTLE_DUP_ADMIT=1 and the ZMW rule, RUSTLE_READ_STRAND unset, layers off.
  ZMW = read name up to and including the ZMW hole number (PacBio `movie/zmw/...`); a coordinate duplicate counts only when its
  ZMW differs from the first record's with that key. "New chains" = spliced chains in the arm absent from base; CAT match =
  exact intron chain of a CAT Liftoff transcript (the r1524 proxy). Base chains kept = r1518 bar 5.
