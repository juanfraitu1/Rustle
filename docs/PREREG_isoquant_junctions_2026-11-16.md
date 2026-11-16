# PREREG r1765: is IsoQuant's junction precision an operating point or a mechanism? (2026-11-16)

User question after row 1764: the same analysis for IsoQuant. Written before any junction-level look at IsoQuant's wrong calls.

Already seen (2026-11-16): whole-set gffcompare intron precision on OR6737, IsoQuant 89.5 vs ours 82.9 (sensitivity 50.2 vs 60.2);
exact chains, IsoQuant 16,677 correct at 0.450, ours cut by reads to that correct count 0.540 (`PREFIX.gtf`, primary pool, today's
defaults). So at chain level IsoQuant's precision is an operating point; the junction level has not been looked at.

Tool `bench/denovo_eval/isoquant_junc.py`, workspace `r1765_isoquant_junc`. Samples: OR6737 first (fast panel: one node, wall limit
1 h), then A119b. Arms: ours `PREFIX.gtf` with `--seed-pool primary` and today's defaults (`artifact_refresh_2026-11-16/{g,h}.pri`);
IsoQuant 4.1.0 de novo `-d pacbio_ccs --polya_trimmed all` (`r1646_isoquant/{g,a}.dn_pt`). IsoQuant does not use secondary
alignments (its log), hence our primary pool. KB3781 is out: no primary-pool run of ours and a different IsoQuant arm there.

- Truth: RefSeq transcript introns (chrom, donor, acceptor), as row 1764.
- Locus unit: the RefSeq gene gffcompare assigns each transcript in its `.tmap` (same gffcompare run for both tools); transcripts
  with no reference gene or class code u / p / x / s are out. Population: genes with a spliced model in both tools.
- Support of a junction: alignment records in the sample BAM carrying exactly that intron, primary and secondary counted separately
  (supplementary skipped); rank by primary + secondary.
- Step B (operating point): at each gene cut our distinct junctions, by support, to IsoQuant's distinct-junction count there (no
  cut when we have fewer). **OPERATING POINT** if our pooled precision at the cut >= IsoQuant's - 0.02 on both samples: stop,
  nothing to model.
- Step C (only on a sample failing B): slip = a wrong junction with an annotated intron in the same tool's set at the same gene,
  both ends within 10 bp, and <= 0.2 of that partner's support; above 0.2 = alternative site; no partner = novel. Splice motif
  recorded. **MECHANISM** if on every sample failing B our slips minus IsoQuant's slips >= half of the wrong-junction gap at the
  matched count. Then count how many of our slips `shiftdrop` already removes before calling it new.
- Anything else: **NEITHER**, report the classes.
- Classes and the motif split of wrong junctions are printed in every case (descriptive).

## Result (2026-11-16): OPERATING POINT on both samples

OR6737:

14,037 RefSeq genes with a spliced model in both tools (`r1765_isoquant_junc/or.txt`, 118 s on one node).

| Set | junctions | correct | precision | wrong | slip | alt site | novel |
|---|---|---|---|---|---|---|---|
| ours, all | 183,104 | 153,090 | 0.836 | 30,014 | | | |
| ours at IsoQuant's count per gene | 150,248 | 139,037 | 0.925 | 11,211 | 192 | 71 | 10,948 |
| IsoQuant | 152,145 | 138,434 | 0.910 | 13,711 | 8 | 8 | 13,695 |

Bar met (+0.016). Descriptive: IsoQuant keeps almost no slips (8 against our 192), and 803 of its wrong junctions are neither
GT-AG nor GC-AG / AT-AC (ours 0). Wrong junctions of both are mostly GT-AG with three or more records, absent from RefSeq.
Chain level, same files: IsoQuant 16,677 correct at 0.450; ours cut to >= 6 reads 17,014 at 0.540 (`final.gtf` >= 5 reads 17,861 at 0.567).

A119b: 20,444 RefSeq genes with a spliced model in both tools (`r1765_isoquant_junc/a1.txt`, 492 s).

| Set | junctions | correct | precision | wrong | slip | alt site | novel |
|---|---|---|---|---|---|---|---|
| ours, all | 302,848 | 203,681 | 0.673 | 99,167 | | | |
| ours at IsoQuant's count per gene | 201,815 | 165,165 | 0.818 | 36,650 | 311 | 222 | 36,117 |
| IsoQuant | 217,700 | 159,154 | 0.731 | 58,546 | 14 | 24 | 58,508 |

Rule: ours at the cut >= IsoQuant - 0.02 on both samples: met (+0.016 OR6737, +0.087 A119b). Step C not required. Descriptive, both
samples: IsoQuant keeps almost no slips (8 / 14 against our 192 / 311), which is under 2% of our wrong junctions at the cut; it keeps
junctions that are neither GT-AG nor GC-AG / AT-AC (803 / 5,751; ours 0) and, on A119b, 17,687 wrong junctions with one or two
records (ours 3,284). Chain level A119b: IsoQuant 18,903 correct at 0.224; ours cut to >= 12 reads 19,499 at 0.414 (`final.gtf` 19,314
at 0.504). Nothing to model from IsoQuant.

Limits: gene unit = gffcompare tmap (transcripts outside RefSeq genes are out); per-gene cut; our slips not checked against
`shiftdrop`; `final.gtf` cut was fixed on these two samples.
