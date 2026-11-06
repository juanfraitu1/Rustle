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
