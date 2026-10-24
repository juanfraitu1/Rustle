# Pre-registration: O3 -> O2 with `--semi --guard` and the site homology check (P14), 2026-10-24

**Status:** committed before any P14 family is simulated or run.

**Method.** `o3port.py --semi --guard` with SITE_FLANK = 5 (register 1396): P13's method, plus a distinguishing site counts only
when both of the read's alignments are gap-free over 5 bases on each side (P13 failed S2 on two reads where an indel shifted one
alignment).

**Panel.** P10 panel families 31-40 by name (never used for O2): GWFAM31, GWFAM312, GWFAM349, GWFAM35, GWFAM363, GWFAM378, GWFAM380, GWFAM387, GWFAM389, GWFAM395.

**Substrate, baseline, decision rules (S1, S2, S3, power, E1) and the reported limit:** identical to P13
(`docs/PREREG_o2_guard_p13_2026-10-24.md`, sections Substrate to Reported; 90 runs). Run
once; nothing changed after looking.
