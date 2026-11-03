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

## Amendment 2 (2026-11-03, before any arm D number; arm D job 4434265 cancelled at 3 min, rerun after this)
The detector's own log on arm A (batched, 212 s) reached 1 of the 69 tandem loci: only 12 have an admissible junction. Reading
CDK11A/B: the crosser `DN_chr1_1084111_5` (3 reads) has a 17-bp first exon INSIDE CDK11A's body (CDK11A's transcripts end
9 kb past the junction donor), so "UP ends before J" is empty; ZNF587/587B likewise. Rule T2's side blocks are redefined: UP(J) =
same-strand non-users of J sharing an exonic base with a crosser's exons upstream of J and NONE downstream; DOWN(J) the mirror.
A block may extend past J. Everything else (MINORITY, HOMOLOGY, the cut, the bars, arms D and E) is unchanged. The first
definition's reach (1 of 69) is recorded as part of the result.
Clarification to Amendment 2 (same hour, before any arm number): with blocks allowed past J, UP and DOWN of an ordinary gene's
isoforms reach into each other's exons and the homology test returns identity 1.000 on shared genomic sequence (9,254 cuts
genome-wide, cov 0.03-0.3 at the tandem loci). UP and DOWN must share no exonic base; otherwise J is not a tandem boundary.

## RESULT (2026-11-03; detector on arm A under Amendment 2 + clarification, `/anvme/workspace/v138da13-rustle/fusedgw/D/{units.tsv,report.tsv,dbg3.txt}`)
Detector wall 145 s (bar 4 PASS after the bucket index; the first run exceeded 30 min, bar 4 FAIL as first recorded). 33,342
admissible junctions, 6,950 pass MINORITY, **42 pass HOMOLOGY** (74 transcripts, 101 cut rows). Reach on the 69 spliced tandem
loci: 28 have an admissible minority junction, **5 get a cut** (GAGE1/12D 0.993, HSPA1A/B 0.986, CDK11A/B 0.968, TBC1D3D/E 0.889,
UGT2B15/17 0.865); the ZNF pairs sit at identity 0.87-0.88 with coverage 0.15-0.36 of the shorter exon union (they share the
KRAB / zinc-finger exons, not the gene). Cut audit of the 101 cut rows: **68 (67.3%) inside one RefSeq coding gene**, 6 between two
genes, 27 with a side on no coding gene.
**❌ FAIL on bars 1 and 3 from the detector's own output** (at most 5 of 69 loci can change; 67.3% > 10%); arms D and E were NOT run
(a node-hour would only confirm bar 1 and measure bar 2 on a list that already fails bar 3). Rule T2 cuts inside single genes
because disjoint isoform blocks of ONE gene (a 5' and a 3' isoform group) are homologous to themselves wherever a gene carries
internally repeated exons; the gate does not distinguish a tandem paralog from a paralogous exon pair.
Closing reading for the user's question: on A119b the near-identical tandem-pair class that any minority + self-homology cut can
reach is 5-7 spliced loci (r1516: 6; here 5), the 4 intronless blobs are O2's tied-read problem, and the 50+ ZNF / OR / PCDHB
"tandem" fusions are readthrough between diverged paralogs that the family edge definition itself does not join. No rule here is
adoptable; the annotation-guided units (r1495 B) remain the only effective cutter for that class.
