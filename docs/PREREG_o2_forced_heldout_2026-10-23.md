# Pre-registration: O2 as forced haplotype-to-copy matches, simulated gorilla held-out (P7), 2026-10-23

**Status:** committed before any P7 family is run with the O2 code. Goal milestone 3 (`ThesisVault/Topics/Goal_O2_O3.md`); the
user said go on 2026-10-23.

**Method.** `o3port.py` as of 791c0bcc (default options). An AS-tied read (in >= 2 copies' piles within 10 AS) is handled as follows:
- it goes to the copy its haplotype is matched to in every maximum allele-slot matching (a forced edge);
- or it goes to the certified new copy it belongs to;
- otherwise it is ABSTAIN.

**Panel.** The fresh gorilla families whose present copies produce tied reads in the tie scan
(`bench/o3_graph/o2tiescan.py`, 10 reads per allele; ties only, no O2 result looked at): GWFAM168, GWFAM224, GWFAM328, GWFAM430.
GWFAM415, 407, 28 and 401 are development (P1) families and were used for development.

**Simulation.** `bench/o3_graph/o2sim.py`, 100 reads per allele, seed 1. Two cases per family:
- control: nothing erased;
- missing: the copy with the most tied reads in the control is masked and left out of the table.

**Truth and baseline.** Truth is each read's origin. The baseline is minimap2's primary.

## Decision rules (pooled over the 8 runs)
- **O2-1:** wrong <= 2% of placed tied reads (a read of a copy absent from the table is right only on a certified new copy).
- **O2-2:** fixes > breaks against minimap2.
- **Verdict.** O2 passes iff O2-1 and O2-2 pass. With nothing placed, O2-2 fails (0 > 0 is false).
- Run once; nothing is changed after looking.

**Development** (P1 families, same simulation): 0 wrong of 24 placed; 1 fix, 0 breaks; 1,196 of 1,220 tied reads abstained
(identical copies with free slots: nothing forced).
