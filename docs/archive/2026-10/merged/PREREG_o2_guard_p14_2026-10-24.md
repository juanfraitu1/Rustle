# Pre-registration: O3 -> O2 with `--semi --guard` and the site homology check (P14), 2026-10-24

**Status:** committed before any P14 family is simulated or run.

**Method.** `o3port.py --semi --guard` with SITE_FLANK = 5 (register 1396): P13's method, plus a distinguishing site counts only
when both of the read's alignments are gap-free over 5 bases on each side (P13 failed S2 on two reads where an indel shifted one
alignment).

**Panel.** P10 panel families 31-40 by name (never used for O2): GWFAM31, GWFAM312, GWFAM349, GWFAM35, GWFAM363, GWFAM378, GWFAM380, GWFAM387, GWFAM389, GWFAM395.

**Substrate, baseline, decision rules (S1, S2, S3, power, E1) and the reported limit:** identical to P13
(`docs/PREREG_o2_guard_p13_2026-10-24.md`, sections Substrate to Reported; 90 runs). Run
once; nothing changed after looking.

## Result (run once, 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2p14`): PASS
All 90 runs completed.

| arm | runs | tied | placed | wrong | MAPQ >= 1 placed / wrong |
|---|---|---|---|---|---|
| control | 30 | 9,553 | 2,496 | 4 | 8,521 / 885 |
| missing | 30 | 5,116 | 935 | 0 | 4,162 / 880 |
| extra, flagged (guarded) | 8 | 1,476 | 0 | 0 | 1,457 / 314 |
| extra, silent (limit) | 22 | 12,556 | 3,988 | 1,585 (39.7%) | 11,458 / 4,965 (43.3%) |

- **S1 PASS** (4 / 3,431 = 0.12%). **S2 PASS** (0 on X'). **S3 PASS** (0.12% < 13.9%). **Power PASS.** **E1 PASS** (0 placed).
- **Verdict: PASSES.**
- **The limit (reported):** O3 flags the extra copy in 8 / 30 `extra` runs; in the 22 silent runs O2 is about as wrong as minimap2
  (39.7% vs 43.3%): an extra copy that fits a free allele slot is invisible to RNA (Theorem A) and misleads both.
