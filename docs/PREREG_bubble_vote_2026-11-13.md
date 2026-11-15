# STEP 0 PER-BUBBLE VOTE ON A FAMILY GRAPH (register 1696, 2026-11-13)

Why: goal "find uses for graphs or variation graphs for the objectives". Rows 1544 and 1649 judged one graph rule only,
path containment (a read is compatible with a path holding all its nodes), and both named the untested follow-up: a
per-bubble vote, O2's certificate in graph form. Row 1649 left 359 reads with "no single path holds the read's nodes".
Written before any vote number was computed; the scorer is first checked on row 1649's own picks (must print 637
assigned, 623 right / 14 wrong, 18 graph-only).

Inputs, all on disk, no new alignment: `/anvme/workspace/v138da13-rustle/r1649_graph/` (8 KB3781 families: seqwish GFA
with paths FAMILY|COPY|KIND|ID, GraphAligner GAF), truth and O2 status `r1642_exonplace/pof_{panel,heldout}.dump.tsv`,
k-mer rule picks `r1647_kallisto/kb.picks.tsv`, `~/runs/gorilla_KB3781/kb.fam.copies.tsv`.

RULE (`bench/o2_contested/bubble_vote.py`), fixed:
- best alignment per read as in row 1649 (longest aligned query span, then identity), >= 0.90 of the read aligned;
- locus = copies with overlapping spans; a locus carries a node when any of its paths holds it;
- informative node = a node of the read's walk carried by some but not all loci of the graph;
- contradictions of a locus = informative nodes of the walk it does not carry;
- assign to the locus with the fewest contradictions iff the runner-up has at least 2 more (margin 2); else abstain.
Margins 1 and 3 are printed for the reader, not judged.

Scoring as row 1649: pick right when its span overlaps a truth copy; decidable molecules only; primary right = the
dump's `al_right` (leftmost-position rule, known to under-count the primary slightly).

Stop rule. PROCEED to a registered held-out row (HG002 on its own haplotype, r1644 assets) only if all three hold:
(a) the vote assigns >= 30 decidable molecules that the k-mer rule AND O2 both leave (row 1649's bar);
(b) wrong rate over all it assigns <= 5% (row 1649's bar);
(c) on the molecules of (a), the vote's wrong rate is lower than the primary's on the same molecules.
Else STOP. Reading registered in advance: if (a) holds and (c) fails, the vote has converged to minimap2's primary
(rows 1642-1644) and adds nothing over it. A STOP says the rule adds nothing over the k-mer rule or the primary on
this panel; it is not evidence against the variation graph as the family model.
Expected: (a) met (ceiling = the 388 decidable molecules the k-mer rule leaves), (c) uncertain.
Limits: one sample, the reference individual, eight families seen by rows 1642-1649; node visits are votes without a
base-quality or CIGAR check; negatives (origin outside the catalog) are reported, not judged.

## RESULT r1696 (2026-11-13, login node, seconds; `/anvme/workspace/v138da13-rustle/r1696_bubblevote/`)
Scorer check on row 1649's picks: 637 assigned, 623 right / 14 wrong, 18 graph-only, k-mer rule 638 / 624: reproduced
(the scorer takes the panel from the families present in the picks: 1,043 decidable molecules against row 1649's
1,026, and 770 at margin 3 where four families assign nothing; assigned, right, wrong and vote-only counts do not depend on it).
Margin 2 (judged): 565 assigned, 553 right / 12 wrong (2.1%; 11 of the wrong in MCL5).
(a) molecules the k-mer rule and O2 both leave: 8 (7 right / 1 wrong); bar 30: NOT met.
(b) wrong rate 2.1%: met. (c) vote 1 wrong of 8, primary 1 wrong of 8: not lower, NOT met.
**STOP.** Not judged: margin 1 gives 765 assigned (741 / 24), 29 vote-only (22 / 7; the primary 24 / 5 on the same
molecules, vote right where the primary is wrong 0, the reverse 2); margin 3 gives 364 (363 / 1), 6 vote-only.
Reading: against containment (row 1649) margin 2 loses MCL215 (59 to 28), MCL5 (41 to 20) and MCL234 (54 to 47)
and gains MCL9 (66 to 87); MCL177 and MCL158 stay near empty under every rule (containment 22 and 17, margin 1 25
and 19 of 136 and 138 reads: ties, not a margin effect). Where the vote adds molecules it does not correct the
primary (0 of 8; 0 of 29 at margin 1, small sets). Containment (row 1649), the
exact k-mer rule (row 1647) and the vote carry the same copy-specific sites of the reference; on this panel no graph
rule of these three shapes over reference-only paths adds reach (eight seen KB3781 families). This is not evidence against the variation graph as the family model:
the graphs here hold one haplotype per copy, so they carry no variation beyond the copies themselves.
