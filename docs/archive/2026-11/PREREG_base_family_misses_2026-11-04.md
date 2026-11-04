# PREREG r1524: recover the base-mode family misses (strand vote, pass-1 drop, catalog-family cut)

User's goal text (2026-11-04), verbatim:
> Mechanism: (1) single-exon models take their strand from the read-orientation vote (RUSTLE_READ_STRAND on) instead of the
> "+" placeholder; (2) find and fix the pass-1 stage that drops a two-read canonical intron chain with no competing model
> (PRAMEF8, chr1:12700000-12740000, two identical MAPQ 22 primaries, polish off, primaries only); (3) the fused-model cut
> from r1523 runs on catalog families (x.fam.copies.tsv), not truth genes, as a final.gtf layer.
> Panel: r1521, 420 primates-level families, A119b, judged under T, plus the 34 missed members with a >= 2-read local chain
> (~/famcmp/assemblable.tsv).
> Bars: of the 34, >= 20 gain a same-strand model; complete families >= 204 + 8; fused families <= 40; ordinary-member hit
> rate >= 0.819 (base value); genome-wide spliced chains identical to base (bar 5 of r1518, >= 0.99); wall: one genome-wide
> assembly run.
> PASS -> (1) and (3) default on in final.gtf layers, (2) fixed in base.
> FAIL -> each part judged alone; a part that passes its bars ships alone.

## Fixed before any number
- Base = fusedgw/A/x.gtf (r1518 arm A, layers off). Arm = one genome-wide assemble run with (1) on and (2) fixed, layers off,
  then (3) applied to its GTF; the cut layer is also scored on base alone so each part is judged alone.
- "Same-strand model" for the 34 = any model (any exon count) on the gene's strand overlapping its exons (r1521 Amendment 2).
- The 34 are listed in ~/famcmp/assemblable.tsv (ours = 0, maxP_chain >= 2); the 3 tied-only members are not in the 34.
- Bar 5: spliced intron chains of the arm GTF vs base, identical share >= 0.99 of base chains (bench/seed_unspliced_score.py bar 5).
- Panel scorer: bench/famcmp/family_panel.py, hit = any model (default), definition T.
- (2) is a code fix only if its cause is found; the fix must keep the 17 family tests and add one test from the PRAMEF8 reads.

## RESULT (2026-11-04; jobs 4435754 combined arm, 4435820 part (1) alone on the pre-fix binary da6c15ee; ~/r1524/panel.out, panel_rs.out)
Combined arm (fix + RUSTLE_READ_STRAND=1, layers off), then the catalog cut:
| bar | base | arm | arm + cut | verdict |
|---|---|---|---|---|
| 1: of the 34, same-strand model | 0 | 27 | 27 | PASS |
| 2: complete families (T) | 204 | 229 | 229 | PASS (>= 212) |
| 3: fused families | 49 | 55 | 53 | FAIL (<= 40) |
| 4: ordinary hit rate | 0.819 | 0.861 | 0.861 | PASS |
| 5: base spliced chains kept | - | 0.9800 | - | FAIL (>= 0.99) |
T-only hit rate 0.690 -> 0.759. Transcripts 236,004 -> 297,472: 68,749 new spliced chains (59,597 at reads = 2) with an exact
CAT chain match of 0.010 vs 0.051 for base's own 2-read chains; 4,376 base chains lost to the shifted polish (ISM collapse and
quantiles). Verdict as registered: FAIL (bars 3, 5). Parts judged alone:
- (1) strand vote alone (job 4435820): bar 5 1.0000 (every spliced chain identical), complete 204 -> 209, ordinary 0.819 -> 0.829,
  fused 49 -> 49, T-only 0.690 -> 0.678 (one member), 3 of the 34 (SLC35G3, CENPVL3, GMCL2), transcripts 236,004 -> 233,011
  (single-exon models 16,501 all '+' -> 13,508 = 10,405 '+' / 3,103 '-'; the 2,993 fewer single-exon models are unexplained
  and the r1518 single-exon bars were not captured for this arm). Passes bars 4, 5; fails 1, 2, 3 as written. A clean, cheap
  gain with no spliced cost; it does not ship under this row's bars.
- Wall: two genome-wide runs (2,462 s and 2,425 s), the second required by the FAIL branch to judge (1) alone; the registered
  bar said one.
- (2) the duplicate rule, by difference (combined minus (1)): +20 complete families and 24 of the 34, bought with 68,749 chains at
  1% CAT match and bar 5 at 0.98. As implemented (any chain reaching the floor through coordinate duplicates) it is too broad:
  most coordinate duplicates are duplicates of 1-read chains. FAIL; the PRAMEF8 class needs a narrower admission (e.g. duplicates
  from different ZMWs only, or the rule confined to catalog-family loci), registered on its own.
- (3) catalog-family cut alone on base: fused 49 -> 47 (the truth-gene cut gave 36), every hit rate unchanged. FAIL (<= 40):
  the catalog's copies (n_reads >= 2) cover far fewer fused pairs than the truth genes do.
Follow-ups: (1) as a registered default-on candidate with its own bars (bar 5 + families + the r1518 single-exon bars);
(2) narrowed; (3) stays opt-in bench.
Housekeeping: pre-fix worktree /anvme/workspace/v138da13-rustle/wt_prefix (binary ~/rustle_target_prev) can be removed.

### After the numbers (housekeeping, no result changed)
Part (2) is opt-in from commit 89bb3fe1: `RUSTLE_DUP_ADMIT=1`; unset, pass 1 is byte-identical to the pre-fix binary. The
default-on version measured above (commit 5674ae49) never reached a shipped run.
