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

## Step 0 RESULT (~/r1536/step0_k1.py logic in step0_k1.out, step0_k1.tsv)
| class | members |
|---|---|
| K1 non-canonical junction in the best own chain | 7 (SPANXA2 53 reads TT-AG; SIGLEC10, RFPL1, OR2L5, FAM86B1, CHST6, RHD) |
| K2 best own chain < 2 reads (every chain unique; the exact-chain floor) | **31** |
| K3/K4 canonical best chain >= 2 reads, no model (window runs pending) | 14 |
K2 >= 10: stop rule not fired. K2 was called "moot" in Amendment 1 because the floor is 2; it is not moot: the floor is met
by no single chain when reads spread over distinct chains (the §6m5 defect; the code's own `--read-isoform-k` widening
admits a chain when every junction has >= k reads). The arm is that knob, measured first on the K2 members' junction support.

### K2 / K3 / K4 measured (~/r1536/k2_junction_support.out, k34_windows.out, k4_duplicates.out)
- K2 (31): 1-6 spliced own reads each, every chain unique. Per-junction widening would admit a canonical chain at k >= 3 for
  5 members, at k >= 2 for 11; 20 have no junction carried by 2 reads.
- K3 (polish): 0 of 14. Windows with `--assembly-polish mono` (ISM / isoform-fraction stages off) change no member.
- K4 (14): every best chain (2-3 reads) is ONE coordinate set from 2-3 different ZMWs: the r1524/r1525 duplicate class
  (RUSTLE_DUP_ADMIT, opt-in only after 68k chains at 1% CAT genome-wide). Closed there; not re-registered.
- One window-only model (PRB2, primaries only, 3 reads) is noted, not pursued.

## Arm (as the step 0 rule obliges: K2 >= 10 -> the widening knob; expected reach <= 5 of 52 at k = 3)
Genome-wide A119b `assemble` with RUSTLE_READ_ISOFORM_K=3, layers off, same binary and AS table as L1; base = fusedgw/L1/x.gtf
(the unlayered GTF of the same binary; gate: L1 x.gtf chains kept >= 0.99 vs R1). Bars (goal row): complete >= 209, fused <= 49,
ordinary hit rate >= 0.829, T-only >= 0.678, catalog pairs precision >= 0.802 and recall >= 0.740 (PAIRSCORE_SAME_SYMBOL=ignore),
new spliced chains with an exact CAT match >= 0.05 (r1531 bar), wall <= 2 h. K2 members gaining a model are reported.
