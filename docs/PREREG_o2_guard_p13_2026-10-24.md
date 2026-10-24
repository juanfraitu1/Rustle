# Pre-registration: O3 -> O2 with `--semi --guard` on fresh gorilla fibroblast families (P13), 2026-10-24

**Status:** committed before any P13 family is simulated or run.

**Method.** `o3port.py --semi --guard` as of 500819a1 (registers 1391, 1394): P12's O2, plus O2 abstains on every tied read
touching a copy where O3 flags a putative copy (or a copy that putative haplotype could take a slot on).

**Development** (GWFAM248, GWFAM415, case `extra`, d 2 / 3): flagged run GWFAM248 d3 91% -> 0 wrong; the silent run GWFAM248 d2
stays 40% wrong (nothing flagged; Theorem A).

**Panel.** P10 panel families 21-30 by name (never used for O2): GWFAM259, 271, 277, 289, 29, 293, 295, 297, 302, 308.

**Substrate.** `o2dup.py`, `O2DUP_CONV=1`, `O2DUP_PORT_ARGS="--semi --guard"`, `O2DUP_CASES=control,missing,extra`, 100 reads per
allele, seed 1, d in {1, 2, 3}: 90 runs. `extra` = control plus a third near copy X* (d substitutions, seed 3) the reference lacks;
X*'s reads are right only on a certified new copy. A failed run is reported and skipped.

**Baseline:** minimap2's primary where MAPQ >= 1, on the same tied reads.

## Decision rules
- **S1:** control + missing pooled: wrong <= 2% of placed. **S2:** 0 tied reads on X' in missing runs.
- **S3:** control + missing pooled: wrong rate < minimap2 MAPQ >= 1's. **Power:** >= 100 placed in control + missing.
- **E1:** `extra` runs where O3 flags a putative copy on X or X' (the guarded runs): O2 wrong <= 2% of placed (0 placed passes).
- **Verdict:** pass iff S1, S2, S3, power and E1. Run once; nothing changed after looking.

**Reported, not decided (the limit):** the number of `extra` runs with no flag (silent), and O2's and minimap2's wrong rates there.

## Result (run once, 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2p13`): FAILS S2
All 90 runs completed.

| arm | runs | tied | placed | wrong | MAPQ >= 1 placed / wrong |
|---|---|---|---|---|---|
| control | 30 | 9,857 | 3,373 | 2 | 8,682 / 891 |
| missing | 30 | 5,301 | 1,252 | 2 | 4,207 / 886 |
| extra, flagged (guarded) | 7 | 1,184 | 0 | 0 | 1,180 / 182 |
| extra, silent (limit) | 23 | 13,188 | 5,654 | 2,254 (39.9%) | 11,947 / 5,154 (43.1%) |

- **S1 PASS** (4 / 4,625 = 0.09%). **S3 PASS** (0.09% < 19.7%). **Power PASS.** **E1 PASS** (guarded runs: 0 placed).
- **S2 FAILS: 2 X|B reads on X' in missing runs** (GWFAM293 d2 r079, GWFAM297 d2 r098). Cause: an indel beside a distinguishing
  site shifts the read's two alignments by 1 bp, so adjacent positions compare non-homologous reference bases and count as sites.
- **The limit (reported):** O3 flags the extra copy in 7 / 30 `extra` runs; in the 23 silent runs O2 is about as wrong as minimap2
  (39.9% vs 43.1%).
- **Verdict: DOES NOT PASS.**
