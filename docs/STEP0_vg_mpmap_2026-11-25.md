# STEP 0: READS MAPPED TO THE FAMILY GRAPH WITH VG'S SPLICE-AWARE MAPPER (register 1945, 2026-11-25)

Why: user, 2026-11-25, after rows 1943 and 1944 (projection through the graph replaces neither the secondary alignments nor
the second assembly): "Try vg on reads directly". Row 1649 mapped reads to per-family graphs with GraphAligner (no splice
model) and found path compatibility = the exact k-mer rule (18 molecules beyond kallisto and O2, bar 30). `vg mpmap` is the
one splice-aware graph mapper; it was built and evaluated on short reads. This step asks whether it runs on Iso-Seq reads
and, if so, whether it assigns reads the other two leave.

## Measurement (KB3781, row 1649's eight families, graphs and reads, unchanged: `r1649_graph/MCL*.{gfa,reads.fa,gaf}`)
- Graph: the family's seqwish GFA (paths = each copy's genomic, transcript and catalogue sequences), nodes chopped to 256 bp
  for the index (`vg mod -X 256`), `vg index` (xg, GCSA k 16 after `vg prune -r`).
- Reads: `vg mpmap -n RNA -l long -A -F GAF`; cap 20 minutes per family on 8 threads.
- Rule, the same for both mappers (`bench/o2_families/graph_pick.py`): best alignment per read, >= 0.90 of the read aligned;
  the read is compatible with every path holding all its nodes; assigned iff those paths belong to ONE copy index.
  (Row 1649 merged copies of one locus through a copies table; this rule does not, for both mappers alike.)
- kallisto's assignments: `r1647_kallisto/kb.picks.tsv`.

## Stop rule (fixed before the run)
PROCEED (score the vg-only reads against the read truth, then a registered row) only if all three hold over the families
built:
- (i) it runs: >= 0.80 of the reads get an alignment over >= 0.90 of their length;
- (ii) reach: >= 30 reads assigned by vg that GraphAligner and kallisto both leave;
- (iii) agreement: of the reads vg and GraphAligner both assign, >= 0.95 on the same copy.
Else STOP. A family over the cap or failing to index is reported as not built and counts against (i).

Prediction (mine): STOP on (i) or (ii). The `long` preset is meant for reads of a few hundred bases; these reads average
about 900-3,000. Where it aligns fully it will repeat GraphAligner's node sets, as GraphAligner repeated the k-mer rule.

Outputs: `/anvme/workspace/v138da13-rustle/r1945_vg_mpmap_step0/`.

## RESULT (2026-11-25, jobs 4504604 and 4504623; `r1945_vg_mpmap_step0/run.4504604.out`, `posthoc.4504623.out`)
Six of eight families built (15-44 s each); MCL5 and MCL9 failed at the GCSA index (`vg prune -r` with the paths restored
filled the node's temporary disk).
| bar | value | verdict |
|---|---|---|
| (i) reads aligned over >= 0.90 of their length | 848 of 855 on the six families; with the 257 reads of the two not built 848 of 1,112 (0.763) | FAILS as registered (bar 0.80) |
| (ii) assigned by vg, left by GraphAligner and kallisto | 27 | FAILS (bar 30) |
| (iii) same copy where both graph mappers assign | 536 of 536 | holds |

**STOP.** Prediction (STOP on (i) or (ii)) held, but my reason for (i) was wrong: `vg mpmap -l long` aligns these reads
end to end where the graph indexes; (i) fails on the two graphs that did not index.
- vg assigns 566 reads, GraphAligner 539, kallisto 537; the 27 extra are 18 in MCL239 (vg aligns reads GraphAligner covers
  under 0.90) and 8 in MCL115. No truth was applied to them (the gate was not passed).
- Post hoc (plain `vg prune`, no bar): MCL9 builds, 66 assigned = GraphAligner's 66 (65 of 65 on the same copy, 0 new).
  MCL5 (6,362 nodes, 107 sequences) indexes but 98 reads do not finish in 20 minutes on 8 threads (no distance index).
- Reading: a splice-aware graph mapper gives the same copy as GraphAligner read for read; the reads no path holds
  (`abstain_no_path` 121 of 138 in MCL158, 113 of 136 in MCL177) stay unassigned under either mapper. As in row 1649, the
  limit is what the paths can tell apart, not the mapper. The largest family is also where the tool does not finish.
