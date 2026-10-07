# The locus-model arms read at two levels of the family hierarchy: MCL clusters (M) and components of the pre-MCL graph (C) (2026-10-07)

**Status: DEV, human CHM13 / CAT-Liftoff ideal windows (NPIP chr16, TBC1D3 chr17, two read replicates each), simulation, circular by construction, the same windows as `docs/LOCUS_UNITS_2026-10-06.md`. Protocol: Amendment 2 of `docs/PREREG_locus_units_2026-10-06.md` (commit 5019106d, written after the level-M results and before any component-level number or F1 arm existed). Tools: `bench/entangled/{hier_clusters.py (2 unit tests), pool_levels.py, run.sh f1|f1q|hier}`. Products: `/mnt/linuxdisk/tmp/entangled_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/{hierC_*,scoreC_*,famscoreC_*,score_f1*}`.**
User direction (2026-10-07): 'lets do 1 and 4' (4 = back to the ideal windows: a hierarchy-aware family criterion, then the evidence-based fusion cut combined with the primary-first representative).

## Answer

1. **The registered level-M verdicts hide real differences.** Reading the same arms at level C (the connected components of the pre-MCL homology graph, a coarsening of every MCL partition: nested by construction) gives pooled IDEAL-FOUND (R plus E) of D 60, P 67, PC 67, **Q 64**, S_Q 62, F1Q 64, S_D 55, F1 60. Every re-run of `mcl_families --dump-graph` reproduced the arm's registered clusters byte for byte (32 of 32 arms-runs valid).
2. **Rule 1 (primary-first representative, arm Q) is a clean repair at level C**: it changes exactly four gene-runs against D and nothing else, NPIPA2 in both NPIP replicates and TBC1D3P4 / TBC1D3P3 (E3 and IDEAL-FOUND 0 -> 1), no loss anywhere; reachable copies found 48 -> 52 of 54, entangled 12 -> 12, cluster precision of K* unchanged (.78), Compara F unchanged. The NPIP replicate-2 'loss' of three copies at level M is a partition split (they are in K*_C), not an evidence loss.
3. **Primaries-only seeding (P, PC) is a net gain at level C and a trade-off per copy**: reachable 53 of 54 (both families YES under the registered rule), entangled 14 of 24 (level M: 7). Of the 10 E4 failures at level M, all 10 are partition splits. The price: cluster precision of K*_C falls to .67 (D .78), and per copy it repairs NPIPA2, NPIPA1, PKD1P6, NPIPA7 (both replicates), NPIPB3 (replicate 1) and TBC1D3P4 / P3, but loses NPIPB4 (both replicates), AC138969.1 and NPIPB13 (replicate 2; E3 breaks: the representative is no longer the copy's chain); NPIPA8 and AC138894.1 change flags without changing status (neither is found in D). Net +7 (R +5, E +2).
4. **F1 (evidence-proven bridge cuts without f1v2's share rule) changes nothing for any target copy at either level**: it cuts more bridges (NPIP 5 and 9 against f1v2's 2 and 1; TBC1D3 2 and 1 against 1 and 1) but zero per-copy flags differ from D in four runs, and NPIPB14P (readthrough as representative), AC138894.1 and PKD1P6 stay failed. F1Q equals Q.
5. **Under the registered component-level rule P, PC and Q are CANDIDATE_C, and all three are LEVEL-SPECIFIC** (at level M none is a candidate), which by the amendment 'earns nothing by itself'. S_D, S_Q, F1 and F1Q are not candidates (F1Q fails the Compara F clause by .012 in NPIP replicate 1).

## Pooled results (four runs; R = 54 reachable gene-runs, E = 24 entangled; R / E found, then R: E2 / E3 / E4)

| arm | level M R / E (R+E) | level C R / E (R+E) | R at C: E2 / E3 / E4 | NPIP R found r1 / r2 (bar 13) at C | TBC1D3 (bar 12) at C | verdict N / T at C | CP of K*_C (min) | Compara F at C NPIP / TBC1D3 |
|---|---|---|---|---|---|---|---|---|
| D | 48 / 12 (60) | 48 / 12 (60) | 50 / 50 / 54 | 12 / 12 | 12 / 12 | NO / YES | .78 | .754 .754 / .775 .805 |
| P | 53 / 7 (60) | 53 / 14 (67) | 53 / 53 / 54 | 14 / 13 | 13 / 13 | YES / YES | .67 | .767 .767 / .775 .805 |
| PC | 53 / 7 (60) | 53 / 14 (67) | as P | 14 / 13 | 13 / 13 | YES / YES | .67 | as P |
| **Q** | 49 / 12 (61) | **52 / 12 (64)** | 52 / **54** / 54 | 13 / 13 | 13 / 13 | **YES / YES** | .78 | .754 .754 / .775 .805 |
| S_D | 38 / 9 (47) | 45 / 10 (55) | 51 / 47 / 54 | 11 / 12 | 11 / 11 | NO / NO | .79 | as D |
| S_Q | 45 / 10 (55) | 52 / 10 (62) | 52 / 54 / 54 | 13 / 13 | 13 / 13 | YES / YES | .78 | as D |
| F1 | 48 / 12 (60) | 48 / 12 (60) | as D | 12 / 12 | 12 / 12 | NO / YES | .78 | .774 .754 / .775 .805 |
| F1Q | 49 / 12 (61) | 52 / 12 (64) | as Q | 13 / 13 | 13 / 13 | YES / YES | .78 | .742 .754 / .775 .805 |

E4 failures at level M (R + E, found-or-not): not in K*_M / of which no edge path to the family (not in K*_C) / partition split (in K*_C): D 2 / 2 / 0; P and PC 10 / 0 / 10; Q 5 / 2 / 3; S_D 13 / 3 / 10; S_Q 11 / 2 / 9; F1 2 / 2 / 0; F1Q 5 / 2 / 3.

## Predictions (Amendment 2)

**H1** (F1 cuts more bridges than f1v2 and finds NPIPB14P in one replicate): the first clause held, the second **failed** (NPIPB14P stays (E2, E3, E4) = (0, 0, 1)). **H2** (F1 loses no reachable copy): held (48). **H3** (F1Q is a CANDIDATE_C): **failed** by the Compara F clause. **H4** (Q finds >= 53 of 54 reachable at level C): **failed**, 52: NPIPA7 in both NPIP replicates stays at locus purity .112 (E2): Rule 1 does not touch the body, whose artifact chains P removes. **H5** (S_D and S_Q are not candidates at level C): held.

## Reading

- Rule 1's effect is the same at both levels (four exactness repairs), but at level M it was masked by one MCL partition flip; the instrument that scores 'the one cluster K*' mixes representation with partition granularity. The decomposition (no edge path against partition split) is the clean way to report E4, and it is parameter-free.
- Structural separators (S_D) lose copies through EVIDENCE (3 of the 13 E4 failures have no edge path at level C, R 45 against 48), not only through splits; with Rule 1 on top (S_Q) the exactness is 54 of 54 but the entangled copies stay at 10.
- Evidence-based bridge cutting (F1 or f1v2) cannot fix the fused entangled copies in these windows: the readthrough PDXDC2P-NPIPB14P carries as many reads as its parents (every simulated transcript has 10 reads) and has no PAS or promoter signature; F1 cuts other bridges and changes no target. Fused-locus copies remain open (NPIPB14P, AC138894.1 with CLN3 and NPIPB7, PKD1P6 with NPIPP1).

## Limits

The same two windows and the same simulator as every earlier arm of this series: nothing here is held out. The component level has no precision guard of its own beyond CP of K*_C; P's CP_C of .67 means a third of its family component is not family. The primary-supported marks of Q come from a second assembly of the same reads. A candidate at level C only earns nothing by the registered amendment; the next registered step is a check on real reads (A119b chr16 / chr17 and gorilla) with the copy-recovery instruments and on ideal windows around other multi-copy regions.
