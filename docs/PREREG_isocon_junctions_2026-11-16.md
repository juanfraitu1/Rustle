# PREREG r1764: is IsoCon's junction precision an operating point or a consensus effect? (2026-11-16)

User question after row 1763: IsoCon scores higher junction precision than ours on the gorilla target panel (KB3781 0.79 vs 0.74;
OR6737 row 1478 0.61 / 0.88 vs 0.58 / 0.83). Is there a junction-specific mechanism to model? Written before any junction-level look.

Tool `bench/denovo_eval/isocon_junc.py`, workspace `r1764_isocon_junc`. Samples: KB3781 (`apetools/GGO_KB3781_fib_mm/ours/ours.gtf`,
today's defaults) and OR6737 (`r1759_cluster2/ours/or.gtf`, assembled before the transfrag-extension default; `isocon_OR6737/isocon.gtf`).
Same target panel (34 targets), copies assigned as in `toolcompare.py`.

- Step A (sanity): reproduce the toolcompare pooled junction precision with its exon-block truth.
- Truth from here on: RefSeq transcript introns (chrom, donor, acceptor), because block boundaries call annotated alternative sites wrong.
- Population: copies with a spliced model in both tools; multi-copy and single-copy targets reported separately and pooled.
- Support of a junction: alignment records in the sample BAM carrying exactly that intron, primary and secondary counted separately
  (supplementary skipped); rank by primary + secondary.
- Step B (operating point): cut our distinct junctions, by support, to IsoCon's distinct-junction count on the same copies (no cut
  when we already have fewer). **OPERATING POINT** if our pooled precision at the cut >= IsoCon's - 0.02 on both samples: stop,
  nothing to model.
- Step C (only on a sample failing B): a wrong junction is a *slip* when an annotated intron at the same copy and in the same tool's
  set has one end identical and the other within 1-10 bp (or both within 10 bp), and the wrong junction has <= 0.2 of that partner's
  support. Splice motif recorded; a canonical GT-AG wrong junction above 0.2 of its partner is an alternative site, not a slip.
  Same classification for IsoCon. **MECHANISM** if on every sample failing B our slips minus IsoCon's slips >= half of the
  wrong-junction gap at the matched count. Then count how many of our slips `shiftdrop` already removes before calling it new.
- Anything else: **NEITHER**, report the classes.

## Result (2026-11-16): OPERATING POINT on both samples

Step A reproduces the block-truth figures (KB3781 0.735 / 0.791). RefSeq intron truth, copies modelled by both tools; the cut is
per copy (our junctions ranked by read support, cut to IsoCon's count at that copy). Outputs `r1764_isocon_junc/{kb,or}.txt`.

| Sample | Class | ours all (n / precision) | ours at IsoCon's count | IsoCon |
|---|---|---|---|---|
| KB3781 | multi (19 copies) | 304 / 0.836 | 265 / 0.906 | 269 / 0.874 |
| KB3781 | single (25) | 406 / 0.877 | 313 / 0.942 | 320 / 0.909 |
| KB3781 | all | 710 / 0.859 | 578 / 0.926 | 589 / 0.893 |
| OR6737 | multi (20) | 349 / 0.785 | 259 / 0.849 | 263 / 0.821 |
| OR6737 | single (27) | 408 / 0.853 | 318 / 0.931 | 318 / 0.925 |
| OR6737 | all | 757 / 0.822 | 577 / 0.894 | 581 / 0.878 |

Rule: ours at the cut >= IsoCon - 0.02 on both samples: met (ahead by 0.033 and 0.016). Step C was not required; the classes were
printed anyway: slips are 4 (KB3781) and 2 (OR6737) of ours against 1 and 2 of IsoCon. Nearly every wrong junction of either tool is
GT-AG with three or more supporting records and no RefSeq intron within 10 bp at the copy: unannotated junctions, not alignment
slips. Uncut, ours has more correct junctions (610 vs 526; 622 vs 510). No junction mechanism to take from IsoCon.

Limits: 34 targets chosen for OR6737 expression; the OR6737 `or.gtf` predates the transfrag-extension default; per-copy cut.
