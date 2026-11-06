# PREREG: chain once on a family graph (row 1552, 2026-11-06 20:40, user-requested parallel lever)

## Question
Read-star chains every molecule against every copy (MCL1: 30,034 molecules x 31 copies, 1,372 s at 16 threads).
A family sequence graph (copies collapsed, divergence as bubbles) lets one graph chaining per read replace 31
chainings, with the copy read off the bubbles the chain passes. Does minigraph's graph chaining give the same
copy calls as read-star at a fraction of the time?

## Set-up (prototype, outside the Rust code; row 1544's inputs on /anvme/workspace/v138da13-rustle/r1544)
MCL1: 31 genomic copy windows (+-2 kb, genomic.fa), 15,828 of its 30,034 molecules recovered (reads.fa);
baseline = gw7 assignments restricted to those reads: 35 assigned, 41 tied, 15,752 ambiguous.
Control: row 1544's `minimap2 -x splice` CLI run on the same reads, 469 s at 16 threads.
Graph: `minigraph -cxggs` over the 31 windows (first = backbone) with `-L` lowered so that PSV-scale divergence
forms bubbles (minigraph's default ignores variants < 50 bp, which would erase the copy signal). Each copy's walk =
its own window aligned to the graph. Reads: `minigraph -cx lr -t 16`. Copy call per read: the copy whose walk
shares the most copy-specific segments with the read's path; assigned iff a unique best with margin >= 1 segment.

## Bars (fixed before any number; loss bound chosen by me, amendable only before the run)
- Time: graph build + read alignment <= 160 s at 16 threads (>= 3x the 469 s control on the same reads).
- Agreement: of the 35 baseline-assigned reads, >= 31 get the same copy (loss <= 4), and <= 1 gets a different
  copy. Of the 41 tied, report how many the graph assigns (ties are the "do not call" class: a high number is a
  false-confidence warning, not a gain).
- Of the 15,752 ambiguous, report the newly assigned count; no bar (unverifiable), reported for the record.
## Prediction
Time passes (minigraph long-read chaining is seconds to a minute on 16k reads). Agreement is the coin flip: if
`-L` cannot be lowered far enough, copies collapse and most reads land on the backbone with no copy-specific
segment (loss near 35 = FAIL, lever closed at the graph-construction step, not the chaining idea). If bubbles
form at PSVs, agreement 28-34 of 35, tie calls 5-15.

## Result arm 1 (2026-11-06 21:20): minigraph collapses the family; closed at construction
r1552 (job 4450720): at -L 50, 10 and 3 the rGFA has ONE segment (31 copies collapsed onto the backbone, no
bubble); build 5-6 s, align 7-8 s, 17,985 GAF lines on 15,828 reads. No copy information exists in that graph, so
the time bar is vacuous and agreement cannot be computed. The predicted construction-step failure.

## Amendment 1 (2026-11-06 21:25, before scoring): arm 2 = seqwish graph + GraphAligner (row 1544 artifacts)
Row 1544's follow-up run (mcl1.f.gfa: seqwish over the 31 units, 86 nodes, 31 P-line walks = the copy walks;
GraphAligner 1.0.20, reads.f.gaf, 15,828 reads aligned, 7,745 end to end) was never scored. Same copy-call rule
(copy-specific segments; unique best with margin >= 1) and the same agreement bars. Time bar: GraphAligner's own
wall on those reads is not in the log; the arm reports agreement only, and a PASS there would send the time
question to a timed rerun before any further claim. This arm uses the UNIT (transcript) sequences, so the
comparison is against the same baseline (gw7) but the graph is of units, not genomic windows.

## Result arm 2 (2026-11-06 21:35): FAIL on agreement; row closed
seqwish graph of the 31 units: 311 segments, 3,271 copy-specific segment memberships. GraphAligner paths scored
by the registered rule: of the 35 baseline-assigned reads 1 same, 27 different copy, 7 uncalled (bar >= 31 same,
<= 1 different: FAIL). Tied reads called 20 of 41; ambiguous reads called 11,767 of 15,752 (false-confidence
warning). Not a label permutation: baseline copy 17 lands on copy18/22/28/8 by read, with best-minus-second
margins of 0-3 segments out of ~65: the graph does not separate copies at ~99% identity (seqwish on a minimap2
asm20 all-vs-all keeps few PSV-scale bubbles). Time bar not evaluated (GraphAligner wall not logged; moot).
Closing: with minigraph the copies collapse to one segment; with seqwish the bubbles exist but do not carry
the copy identity read-star finds in the base-level alignment. A rarity-weighted scoring rule was not tried
(post hoc). Lever closed; the per-base PSV evidence read-star uses has no cheap graph proxy in these tools.
