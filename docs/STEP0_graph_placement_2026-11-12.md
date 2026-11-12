# STEP 0 GRAPH-ALIGNMENT PLACEMENT (register 1649, 2026-11-12)

Why: user request; third option of the 2026-11-12 table and the thesis's model (variation graph). Question for a step 0
only: does aligning reads to a per-family sequence graph place molecules that the exact k-mer rule (r1647) and O2 both
leave, at a comparable wrong rate? Its only edge over k-mers is tolerance to errors and alleles.

Build (per family): sequences = r1643's per-copy set (cat, tx, pre). GraphAligner 1.0.20 has no splice mode, so introns
exist only as transcript paths. All-vs-all `minimap2 -c -X -x asm20` (no `N` in CIGARs), `seqwish -k 31` (minimum
match 31, against intronic repeat collapse), paths = input sequences. Reads: the family's molecules (decidable +
negative) from the r1642 set. `GraphAligner -x vg`, best alignment per read with >= 0.90 of the read aligned.
RULE (`bench/o2_contested/graph_place.py`): a read is compatible with every sequence whose path contains all nodes of
its alignment; assign iff the compatible sequences are copies of one locus.

Families, fixed before the run: KB3781 families with >= 40 decidable molecules in r1642's set whose family holds <= 30
copies, in the order of r1642's table, at most 8; MCL5 (where the primary errs) is in by this rule.
Wall cap: 20 min per family build + align on 8 threads; a family over the cap is reported as not built.
Scorer: locus overlap for the pick, record overlap for the primary.

Stop rule: PROCEED to a registered row only if, over the families built, (a) the graph assigns >= 30 decidable
molecules that kallisto AND O2 both leave, and (b) its wrong rate on all it assigns is <= 5%. Else STOP.
Not judged: negatives (a per-family graph has no outside rival; reported).
Expected: STOP: what it adds over k-mers is reads with errors at non-PSV positions, few on HiFi.

Outputs: `/anvme/workspace/v138da13-rustle/r1649_graph/`.

## RESULT r1649 (2026-11-12, login node, <= 2 s per family; `r1649_graph/{build.log,step0.txt}`)
Families by the rule: MCL115, MCL9, MCL177, MCL5, MCL239, MCL158, MCL215, MCL234 (all built; graphs 23 to 6,362 nodes).
1,026 decidable molecules. Graph assigns 637: 623 right / 14 wrong (2.2%; all 14 in MCL5). kallisto on the same
molecules assigns 638: 624 / 14. Both assign 613, always the same locus.
(a) graph-only molecules (kallisto and O2 both abstain): 18 (15 right / 3 wrong), bar 30: not met. (b) wrong rate 2.2%, met.
**STOP.** The graph alignment reproduces the k-mer rule read for read and adds 18 molecules; kallisto assigns 25 the
graph leaves. Abstentions are mostly "no single path holds the read's nodes" (359 reads over the 8 families).
Negatives: 9 of 86 assigned (not judged).
Limit of this step 0: node-set containment, per-family graphs, one sample. It says nothing against the variation graph
as the model; it says read placement by path compatibility carries the same information as exact long k-mers here.
