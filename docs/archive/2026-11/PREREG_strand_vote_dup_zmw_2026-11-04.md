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

## RESULT part (1): strand vote (arm R1 = job 4435820; ~/r1525/sus_R1.out, ~/r1525/mono_drop.out, ~/r1524/panel_rs.out)
| bar | base | arm | verdict |
|---|---|---|---|
| spliced chains identical | - | 219,503 / 219,503 = 1.0000 | PASS |
| complete families (T) | 204 | 209 | PASS |
| ordinary-member hit rate | 0.819 | 0.829 | PASS |
| CAT single-exon genes hit (r1518 bar 4 as registered) | 1,465 | 1,490 (ratio 1.017) | PASS |
| single-exon drop explained | 16,501 | 13,508 | PASS, see below |
The 3,088 base single-exon spans absent from the arm (all reads >= 6): 3,031 lie inside a same-strand spliced model of the arm
(the mono-shadow polish now reaches models whose '+' placeholder had put them on the wrong strand), 49 absent, 8 moved <= 50 bp.
95 single-exon models gained, all '-'. T-only hit rate 0.690 -> 0.678 (one member), fused 49 -> 49.
PASS -> shipped: `tools/rustle_pipeline.sh` assemble sets RUSTLE_READ_STRAND=1 unless the caller sets it (0 = placeholder).

## RESULT part (2): different-ZMW duplicate admission (arm R2 = job 4435861, 3,284 s; ~/r1525/r2_chains.out)
| bar | value | verdict |
|---|---|---|
| new chains' exact CAT match | 676 / 68,749 = 0.010 | FAIL (>= 0.05) |
| base chains kept | 215,127 / 219,503 = 0.9801 | FAIL (>= 0.99) |
| of the 34 gain a model | 26 | PASS (>= 15) |
Transcripts 236,004 -> 300,583. The arm is the r1524 blanket arm to the chain: the ZMW restriction removes nothing, because
coordinate duplicates from the SAME ZMW do not exist in this BAM (chr20 spliced primaries: 13,519 coordinate duplicates, 0 from the
same ZMW). Identical (start, end, chain) across distinct ZMWs is therefore the duplicate signature itself, PCR duplicates of
one cDNA, and those chains match CAT at 1% against 5% for base's own 2-read chains. The PRAMEF8 class (26 of 34) is real but
inseparable from the duplicates by read identity alone. The switch acts (region run: PRAMEF8 emitted with RUSTLE_DUP_ADMIT=1,
absent with 0). Verdict: FAIL; stays opt-in. A separating signal would have to come from sequence (two duplicates with
different error profiles are two molecules; identical sequences are one), not from names.
Row verdict: (1) PASS and shipped; (2) FAIL, opt-in. The fused-family cut is closed as de novo (r1524).

### After the numbers (housekeeping, no result changed)
- The ship commit b8d29f0d was inert for the default (good) seed pool: the pool branches reassign `seed_env` wholesale. Fixed by
  appending RUSTLE_READ_STRAND after the branches; verified by a pipeline assemble of chrM+chr21 with the default pool
  (single-exon models 401 '+' / 58 '-'; all '+' before).
- The ZMW machinery (name hash per record, seen-map, counter) is removed since it is vacuous here; RUSTLE_DUP_ADMIT stays as the
  opt-in blanket rule, read once per region. Commit 03ccd245 holds the exact code that produced arm R2.
- Walls: A 2,462 s (r1524 combined arm), R1 2,425 s, R2 3,284 s; the ZMW map cannot plausibly explain the 800 s, left unattributed.
- Sequence lever pre-check (2026-11-04, not registered): coordinate-duplicate pairs vs independent same-chain pairs (different
  ends), compared over a shared 1.5 kb window by aligned bases + indels, chr20, 400 pairs each: duplicates identical 391/400,
  independent 397/400; at the 34 misses 146 of 154 duplicate pairs identical. HiFi reads are error-free at this resolution and
  Iso-Seq molecules of one isoform share TSS/TES, so sequence carries no PCR-duplicate signal (the earlier "0 mismatches" came
  from NM being all indels). The sequence rule is closed without a row; RUSTLE_DUP_ADMIT stays the opt-in blanket rule.
