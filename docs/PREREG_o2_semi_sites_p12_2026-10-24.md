# Pre-registration: O2 `--semi` with supported distinguishing sites on fresh gorilla fibroblast families (P12), 2026-10-24

**Status:** committed before any P12 family is simulated or run. Goal milestone 3, user rule "only exact ties are hopeless".

**Method.** `o3port.py --semi` as of 3b11426e (register 1391). A tied read nothing forces: HOPELESS on an exact tie (same AS and
edit distance on its two best copies); otherwise placed on the leader only when it leads by >= 2 distinguishing sites (reference
bases of the two copies differ), each counted only if the read's base is carried by >= 6 reads of the leader's pile at that column.

**Development** (GWFAM248, GWFAM415; d 0/1/2, conversion): 582 placed, 1 wrong (0.17%), 0 on the lacked copy.

**Panel.** P10 panel families 11-20 by name (never used for O2): GWFAM176, 178, 189, 193, 194, 226, 233, 255, 257, 258.

**Substrate, truth, baselines, decision rules:** identical to P11 (`docs/PREREG_o2_semi_p11_2026-10-23.md`): o2dup.py,
`O2DUP_CONV=1`, `O2DUP_PORT_ARGS=--semi`, 100 reads per allele, seed 1, d in {1, 2, 3}, cases control / missing, 60 runs.
- **S1:** wrong <= 2% of placed tied reads. **S2:** 0 tied reads on X' in missing runs. **S3:** wrong rate < minimap2 MAPQ >= 1's.
- **Power:** >= 100 placed. **Verdict:** pass iff all four. Run once; nothing changed after looking.
