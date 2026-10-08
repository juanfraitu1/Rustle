# The current default at NPIP, re-scored at HEAD (2026-10-06)

Protocol and rules: `docs/PREREG_default_rescore_npip_2026-10-06.md` (committed 7c26e2c0, amended 1df2e20b before any DEF or PRE product existed). Runner:
`bench/default_rescore/run.sh`. Products: `/mnt/linuxdisk/tmp/rescore_2026-10-06/` (`verdict.json`, `support.copies.tsv`, `fs_DEF.json`, `fs_PRE.json`, `provenance.txt`).
Binaries: HEAD 8fbaba79 build in `rustle_target_m2/release` (`copy_assign` 87824d91, `mcl_families` a6308244, `family_score` 542923fd; `copy_support.py` d0ea1bfb).
Human A119b chr16, the 25 CAT/Liftoff NPIP copies. **Everything here is DEV** (the development block of NPIP and of every fusion rule): a measurement, not a validation.

## Verdict: PASS (R1 to R4 hold; G0, G1, G3 pass)

| rule | value (DEF = the default) | bar | |
|---|---|---|---|
| R1 strict FOUND, Amendment E, within own nodes | **8 of 25** | >= 6 | ok |
| R2 locus-level E-found | **21 of 25** | >= 18 | ok |
| R3 U2 bipartite F | **0.645** (sens .588, prec .714; pairs 117 of 240 truth, 199 predicted) | >= 0.645 | ok, at the bar |
| R4 own node of NPIPB2 and NPIPB6 | both have one (PRE: neither) | both | ok |

R3 sits exactly on its bar because the HEAD products are byte-identical to the 2026-09-30 products that set it (below), not because of a margin.

Gates: G0 the HEAD scorer reproduces the registered Amendment E numbers on the frozen arms exactly (P 23 / 10 / 19 / 10, GOOD 21 / 6 / 18 / 6, ALL 24 / 2 / 16 / 2: own node, E-found in own nodes, locus level, E-found chr16-wide);
G1 the own-node rule reproduces 50 of 50 stored flags (P and GOOD; the frozen ALL arm has two loci on one span, see the amendment); G3 the HEAD `family_score` reproduces the stored default-arm U2 (.588 / .714 / .645) and Compara (.500 / 1.000 / .667).
G2 (diagnostic): the PRE control, rebuilt by the HEAD binaries, reproduces the registered GOOD arm exactly (21 / 6 / 18 / 6). The scorer is deterministic under two `PYTHONHASHSEED` values.

## DEF against PRE (descriptive; same binaries, same data)

| | PRE (`--bridge-regroup off`, cov 0) | DEF (default: f1v2, cov .70) |
|---|---|---|
| copies with an own node | 21 | 24 |
| strict E-found in own nodes (of 25) | 6 | **8** |
| locus-level E-found (of 25) | 18 | **21** |
| U2 sens / prec / F | .500 / .680 / .576 | .588 / .714 / **.645** |
| Compara chr16 sens / prec / F (17 families) | .444 / 1.000 / .615 | .500 / 1.000 / **.667** |
| Compara CF153 (the NPIP family): genes hit, F | 10 of 19, .690 | 13 of 19, **.812** |
| Soto sens / prec / F | .465 / .917 / .617 | .507 / .923 / .655 |

Copies that change between PRE and DEF: **NPIPB2** and **NPIPB6** go from no node to an own node and are found (their holders were fused with GSPT1 and EIF3CL before the bridge regroup); **NPIPA6** gains an own node and is not found.
The E-found copies in DEF: NPIPB2, NPIPA2, NPIPA1, NPIPA5, NPIPB6, NPIPB8, NPIPB10P, NPIPB11.

## What remains short

- 22 of 25 copies are spliced-expressed under E (PKD1P6-NPIPP1, NPIPB12 and LOC124907808 are not). **14 of those 22 are not found** by the default: the locus representative does not carry the capped start and its first three introns, although a transcript of the same locus does at 13 of them (NPIPB5 at none). The loss is the representative, as `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` found; f1v2 does not repair it (+2 copies, both from node separation).
- Found is 32 % of the copies in own nodes (40 % in the primaries-only arm P before f1v2: P was not re-run on HEAD, so no DEF-versus-P statement is made).
- Representative rules that would act on this (R_J, most junctions) were tested and are opt-in (`docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md`); they were not run here.

## Provenance: HEAD reproduces the 2026-09-30 products byte for byte

`cmp` of the HEAD products against the e163d955 container-headroom products: **identical** for DEF (assembled GTF, `families.gtf`, clusters, loci GFF3, copy table) and for PRE (assembled GTF, clusters, loci GFF3).
So the 2026-10-06 consolidation of `src/` (`2a83a747`) and the later opt-in code changed nothing in the default's chr16 outputs (assemble 24 s / 0.75 GB; families 65 s). This is one input (A119b chr16), not the whole-genome run.

## Not shown

Any held-out NPIP substrate, any gorilla substrate (strict E needs the cap signal, which gorilla lacks), O2, O3. Liftoff recall was not recomputed.
