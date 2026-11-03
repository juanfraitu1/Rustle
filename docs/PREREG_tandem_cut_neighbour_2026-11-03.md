# PREREG TANDEM CUT, NEIGHBOUR BLOCKS (register 1517, 2026-11-03)

**Written before any number of the arm.** Follows r1516 (`docs/PREREG_tandem_cut_2026-11-03.md`): Rule T found a second block
inside the crosser's own gene_id at only 16 of the 69 spliced tandem loci. In 34 of the 63 two-gene loci the second paralog's own
transcripts lie in ANOTHER gene_id (the crosser shares no junction with them) and in 21 the locus is the crosser alone. User
(2026-11-03): "should we preregister that row?" -> yes.

## Rule T2 (fixed now; `bench/tandem_cut.py --mode neighbour`)
Per transcript t with a `reads` attribute, per intron J = [s, e] of t (1-based closed), on t's strand:
1. CROSSERS(J) = every transcript using J (any gene_id). UP(J) = same-strand transcripts NOT using J whose last exon ends < s and
   that share >= 1 exonic base with a crosser's exons upstream of J; DOWN(J) = same-strand transcripts not using J whose first
   exon starts > e and that share >= 1 exonic base with a crosser's exons downstream of J. Any gene_id. J is admissible when UP
   and DOWN are non-empty and every transcript involved has a `reads` value.
2. MINORITY(J): reads(CROSSERS) < min(reads(UP), reads(DOWN)).
3. HOMOLOGY(J): exon-union of UP vs exon-union of DOWN, `minimap2 -c -x asm20`, identity >= 0.80 and aligned span >= 0.50 of the
   shorter (r1516's gate, unchanged).
4. Every crosser is cut at J (row `tid chrom strand junctions label gene`, label `tandem2`). A transcript cut at several J gets the
   union (the list parser merges rows). No recursion is needed: every J is judged on its own. Rule T's cuts are NOT included (T2
   subsumes T wherever T's blocks exist, with the same thresholds).

## Arms (A119b genome-wide, r1495 protocol, `~/fused/cgw.sbatch` with the T2 list as arm D)
A as before. D: `RUSTLE_BRIDGE_REGROUP=f1units RUSTLE_BRIDGE_UNITS_LIST=<tandem_cut.py --mode neighbour on A/x.gtf>`.
Metrics as r1516: merge_audit fused loci, the spliced tandem class (69, the `tandem_fused.py` classification on arm A's loci,
re-identified in D by the RefSeq gene pair), pairscore human coding precision / recall, cuts inside one RefSeq coding gene,
PKD1P6-NPIPP1 present, detector wall. Reported beside: reach (how many of the 69 get a cut) by r1516 class (one-sided /
span-only / two-sided).

## Decision rule (fixed now)
1. Spliced tandem-class fused loci: 69 -> <= 34.
2. Pair precision >= 0.7158 and recall >= 0.7307.
3. Cuts inside one RefSeq coding gene <= 10% of all cuts.
4. Detector wall <= 30 min.
PASS all four -> `--mode neighbour` is the opt-in de novo tandem detector (default unchanged; OR6737 replication before any default
decision). Any FAIL -> recorded; a one-sided cut whose RIGHT block is another gene_id's transcripts stays documented as the
mechanism either way.

## Amendment 1 (2026-11-03, before arm D ran; from r1516's arm C post hoc)
`--bridge-units-list` replaces F1's own bridge nomination (arm C: fused 547 -> 667 with only 19 tandem cuts). Arm D's list is
therefore the UNION of (i) F1's bridge junctions on arm A (`A/x.bridge_junctions.tsv` rows with `bridge` True, every TJ transcript
cut at that junction: what plain `f1units` nominates) and (ii) Rule T2's cuts. A control arm E = list (i) alone (plain f1units)
runs in the same job and is reported beside D so the tandem effect is D minus E; the bars stay as written (vs A). Bar 4 is
already FAILED: the T2 detector on the whole assembly exceeded 30 min (it was still running at 1 h); arm D runs regardless for
bars 1-3.
