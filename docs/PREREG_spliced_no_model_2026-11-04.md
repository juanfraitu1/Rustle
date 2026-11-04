# PREREG r1536: why pass 1 emits no spliced model at the 52 B3 members (step 0 = diagnosis; 2026-11-04)

From r1534: 52 panel members have >= 2 own spliced primaries and no spliced model (SPANXA2 58 own reads, 53 on one intron chain,
MAPQ 13-22, zero models in a primaries-only 400 kb window, ~/r1534/s_prim). Goal row (docs/GOAL_family_assembly_2026-11-04.md).

## Step 0 rule (fixed before any trace)
Per B3 member, in a primaries-only +-200 kb window (`samtools view -F 0x900`, `assemble --contig`, same binary, as-table from the
full BAM), the member is classed by the FIRST knob whose change makes a spliced model on its strand overlap its exons:
K1 junction strictness (RUSTLE_JUNCTION_MAJORITY=1), K2 the chain floor (pass-1 min reads; name read from the code before the
run), K3 polish off (ISM / isoform-fraction drops), K4 none of these (unexplained). Knobs are applied one at a time, never
combined. Stop rule: if no single knob explains >= 10 of 52, the class is closed as "no single rule", no arm. If one does, the
arm is that knob's rule re-registered as a candidate default with the r1521 panel bars (goal row) and the r1530 catalog score.
Wall: windows run in seconds each; whole step 0 < 30 min.

### Amendment 1 (before the 52 are classed; after one SPANXA2 window)
RUSTLE_JUNCTION_MAJORITY=1 does not switch off pass 1's canonical-junction gate (copy_assign_bin.rs:329: "pass-1 already
enforces canonical GT-AG"; majority only tolerates a non-canonical minority next to a canonical majority). SPANXA2's single
junction (139898565-139899212, exactly the CAT exon boundary) reads TT..AG on CHM13, so the gate rejects every one of its 53
reads. K1 is therefore measured directly, without a run: K1 = the member's best own chain (>= 2 reads) carries a junction whose
motif on the member's strand is not GT-AG / GC-AG / AT-AC. K2 is moot (PASS1_MIN_READS = 2). K3 (polish) and K4 are classed
by a window run only for members not in K1. Stop rule unchanged (>= 10 of 52 for one class).
