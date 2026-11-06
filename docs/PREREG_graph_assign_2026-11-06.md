# PREREG: read-to-graph assignment pre-check on MCL1 (row 1544)

Written 2026-11-06 05:45 before any number.

## Question
Is aligning each read once to a family variation graph cheaper than read-star's minimap2 against every copy,
and does it reproduce the copy calls? MCL1 (A119b, chr13 13.98-15.18 Mb, 31 de novo copies, 30,035 reads) is
the worst family: read-star 741 s at 16 threads (gw7), 48 assigned, 64 tied, 29,923 ambiguous.

## Arms (all on one node, 16 threads, 1 h wall, outputs /anvme/workspace/v138da13-rustle/r1544)
1. **CLI control**: `minimap2 -c --eqx -x <preset> --secondary=yes -p 0.3 -N 31 -t 16` exactly as read-star
   runs it, genomic form (copy genomic windows, `-x splice`) and unit form (31 unit sequences, `map-hifi`),
   on the same 30,035 reads. Splits wrapper cost from aligner cost. If the CLI is minutes, the lever is the
   wrapper (identity-preserving) and beats the graph for the speed goal.
2. **Graph**: `minimap2 -cx asm20 -X` all-vs-all on the 31 unit sequences -> `seqwish` -> GFA with one P path
   per copy (minigraph not used: its minimum variant length drops SNP-scale PSVs). `GraphAligner -x vg -t 16`
   -> GAF. Copy call: the set of copy paths containing the read's node path as a contiguous oriented
   subpath; exactly one = assigned, several = ambiguous, none = unaligned.
3. **Reads**: primaries (`-F 0x904`) over the hull, joined to the 30,035 names from gw7's assignments table;
   recovery count reported.

## Bars
- GraphAligner wall (16 threads) < CLI control wall AND < 741 s.
- Every one of the 48 reads read-star ASSIGNED lands on the same copy; zero contradictions.
- Reported, not barred: how many of the 64 tied and 29,923 ambiguous the graph calls uniquely (path
  containment is more liberal than read-star's certificate, so it will likely call more).
Unit form only: GraphAligner is not splice-aware. FAIL on either bar closes the pre-check; PASS opens a row
on the assign stage, which changes outputs and gets its own precision bars.

## Result (2026-11-06 06:05)

Reads: 15,828 of 30,035 recovered as primaries over the hull; the other half have their primary elsewhere
(they enter read-star through secondaries). All timings below are on the 15,828, 16 threads, one node.

| arm | wall |
|---|---|
| CLI control, genomic form (`-x splice`, 31 loci + 2 kb flanks) | 469 s |
| CLI control, unit form (`map-hifi`, 31 units) | 39 s |
| graph build (asm5 all-vs-all, identity >= 0.95, length >= 1 kb; seqwish -k 79): 311 nodes, 11.5 kb for 82 kb of units | 8 s |
| GraphAligner `-x vg` | 14 s |

Speed bar PASS: 14 s vs 39 s (unit) / 469 s (genomic) / 741 s (read-star on 30k). The CLI reproduces
read-star's cost, so the wrapper is NOT the cost: the genomic `-x splice` form is 12x the unit form and is
~90% of read-star's time on this family. First attempts failed and are recorded: seqwish -k 0 on asm20 gave
86 single-base nodes (GraphAligner: 0 seeds); unfiltered asm20 PAF (6,608 lines, identity down to 0.54 from
internal repeats) gave an 825 bp graph.
Call bar FAIL as registered: with the pre-registered rule (read path contained in exactly one copy path), 0 of
the 48 read-star-assigned reads are called (35 "none": the alignment path matches no single copy path exactly,
13 not recovered); 130 of 15,605 aligned ambiguous reads become unique, 17 multi. Contiguous containment is
the wrong rule for HiFi reads with ~0.1-1% errors through SNP-scale bubbles; a per-bubble vote (which copies
agree with the read at each PSV bubble, read-star's certificate in graph form) is the rule a real row needs.
Verdict row 1544: speed PASS, call FAIL as registered; not adopted. Two levers it exposed, cheaper than a graph:
(a) the genomic splice form is 90% of read-star's cost; whether it is needed on families whose units already
explain the reads is a registered question (§6fp two-form rule kept it for NPIP unit 17); (b) half of a
hotspot family's reads are not even primaries there.
