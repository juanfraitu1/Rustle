# PREREG TANDEM CUT (register 1516, 2026-11-03)

**Written before any number of the arm.** User (2026-11-03): "focus on tandem repeats and readthroughs that merge them and try
to develop rules to see if they should be separated" -> "lets pre-register it".

## What is known (descriptive, `bench/annot_compare/tandem_fused.py`, outputs `~/fused/tandem_fused.tsv`, `~/fused/tandem_junction_molecules.tsv`)
r1495 arm A (A119b genome-wide, f1v2 default): 547 loci of `x.families.gtf` join >= 2 coding RefSeq genes that share no exon
(merge_audit). 73 join paralogs (20 Compara truth pairs at the primates level) or same-stem names (53). Three mechanisms:
- 4 single-exon coverage blobs over intronless tandem arrays (USP17L 41 genes / 265 kb, TAF11L 31, OR10A2/4, SPINK7/9); no
  junction exists and the reads are ms-tied (USP17L 67 of 130 primaries) -> O2's problem, OUT OF SCOPE here.
- 69 spliced loci: F1's structural test fails on the best gene-boundary-crossing junction in all but 9 (no private side 35,
  another crosser straddles 16, no crossing intron 9) and the 9 fail its PAS / promoter proofs. The crossing-junction reads
  are real placements: 805 primaries, 8 MAPQ0, 14 ms-tied within 10 (RBMY and TBC1D3 are the tied exceptions). In 34 of 63
  two-gene loci the second paralog's own transcripts lie in ANOTHER gene_id (the crosser only extends the first gene's extent).
- Register metric: those 69 loci hold ~12 expressed truth pairs; the arrays' 1,325 pairs are not in the expressed universe
  (one CAT gene per copy). Refuted before: blunt readthrough drop (loses PKD1P6-NPIPP1), per-transcript relative depth (r808),
  intron cap (r749), read-through filter on blind nodes.

## Rule T (fixed now; `bench/tandem_cut.py PREFIX.gtf GENOME.fa -> units.tsv`, a `--bridge-units-list` detector)
Per input gene_id g of the f1v2 `PREFIX.gtf` (arm A's `x.gtf`), per strand, over the transcripts with a `reads` attribute:
1. Candidate cut points p = the midpoints of every intron of every transcript of g. For p: LEFT = transcripts whose last exon
   ends < p, RIGHT = first exon starts > p, CROSSERS = the rest. p is admissible when LEFT and RIGHT are non-empty, every
   crosser has an intron containing p (a single-exon crosser blocks p), and no crosser lacks a `reads` value.
2. MINORITY(p): reads(CROSSERS) < min(reads(LEFT), reads(RIGHT)), reads = the `reads` attribute summed (F1v2's rule, pooled
   over all crossers, no PAS / promoter proof, straddling allowed).
3. HOMOLOGY(p): the exon-union sequence of LEFT vs that of RIGHT (transcription orientation), `minimap2 -c -x asm20` best hit
   by matches: identity (matches / alignment block length) >= 0.80 and aligned query span >= 0.50 of the SHORTER sequence
   (the family edge definition, `docs/seeded_family_definition.md`). This is the tandem gate: a readthrough between unrelated
   genes never passes it.
4. Choose, among admissible p passing 2 and 3, the one with the largest min(reads(LEFT), reads(RIGHT)) (ties: leftmost);
   record the cut: every crosser is cut at its intron containing p (one `tid chrom strand junctions label gene` row, label
   `tandem:<identity>:<share>`); then recurse on LEFT and on RIGHT as separate transcript sets (depth <= 3). No transcript is
   deleted: f1units turns each crosser into units. Thresholds are fixed here and not tuned.

## Arms (A119b genome-wide, the r1495 protocol, `~/fused/fgw.sbatch` arm function)
- A: r1495 arm A, already run (`/anvme/workspace/v138da13-rustle/fusedgw/A`): fused 547, pairs precision 0.7258 / recall 0.7407.
- C: `RUSTLE_BRIDGE_REGROUP=f1units RUSTLE_BRIDGE_UNITS_LIST=<tandem_cut.py on A/x.gtf>`, assemble + families as A.
Metrics on C's `x.families.gtf`: merge_audit fused loci, the tandem class (paralog + same-stem, the `tandem_fused.py`
classification, spliced only: the 4 blobs are excluded from numerator and denominator), pairscore human coding precision /
recall, and the cut audit: of all cuts, those whose two flanking units overlap exons of the SAME RefSeq coding gene and of no
other (a cut inside one gene). Reported beside: how many of the 69 spliced tandem loci pass HOMOLOGY at all (the gate's reach),
PKD1P6-NPIPP1 present as a transcript in C's families input, wall of the detector.

## Decision rule (fixed now)
1. Spliced tandem-class fused loci: 69 -> <= 34 (a fall >= 50%).
2. Pair precision >= 0.7158 and recall >= 0.7307 (each within 0.01 of A).
3. Cuts inside one RefSeq coding gene <= 10% of all cuts.
4. Detector wall <= 30 min on one node.
PASS all four -> Rule T is an opt-in de novo detector (`tandem_cut.py` -> `--bridge-units-list`), documented beside r1495's
annotation-guided mode; the default stays f1v2 (a default change is the user's decision after OR6737 replication). Any FAIL
-> recorded in the register; the descriptive decomposition stands either way.
