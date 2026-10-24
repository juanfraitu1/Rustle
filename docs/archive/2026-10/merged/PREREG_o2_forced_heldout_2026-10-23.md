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

## Result (run once, 2026-10-23)
| run | tied | placed | wrong | abstain | fixes | breaks | minimap2 wrong |
|---|---|---|---|---|---|---|---|
| GWFAM168 control | 200 | 0 | 0 | 200 | 0 | 0 | 200 |
| GWFAM168 missing | 0 | | | | | | |
| GWFAM224 control | 200 | 0 | 0 | 200 | 0 | 0 | 200 |
| GWFAM224 missing | 0 | | | | | | |
| GWFAM430 control | 2 | 2 | 0 | 0 | 2 | 0 | 2 |
| GWFAM430 missing | 12 | 9 | 0 | 3 | 9 | 0 | 12 |
| GWFAM328 control (27 copies) | 42 | 42 | **24** | 0 | 15 | 0 | 39 |
| GWFAM328 missing | not run: exceeded the heavy-job time limit twice | | | | | | |

**Pooled over 7 runs:** 53 placed, 24 wrong (45%); 26 fixes, 0 breaks.
- **O2-1 FAILS** (bar <= 2%). The missing 8th run could only change this by placing more than 1,100 tied reads with no new error.
- **O2-2 passes.**
- **Verdict: O2 as forced matches DOES NOT PASS.**

**What it shows:**
- In two- and three-copy families the rule is safe but mostly abstains: identical copies with free slots force nothing.
- In the 27-copy family, forced edges were wrong 24 of 42 times. With many near-identical copies, haplotypes merged across piles
  (the read-overlap merge) can be forced onto the wrong copy.
