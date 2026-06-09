# VG mode regresses below baseline — the leak that cancels VG's wins

**Date:** 2026-06-09. Surfaced while characterizing the no_overlap structural-lever candidates.

## Finding

Genome-wide FSM (annotated-isoform recovery), baseline `-L` (strand flag OFF) vs the VG-run
(`--vg --vg-snp` + RUSTLE_VG_TANDEM + RUSTLE_VG_DECISIVE_GATE, the genome-wide protocol config):

| | FSM |
|---|---:|
| baseline total | 24,374 |
| VG total | 24,373 |
| **VG regressions** (baseline FSM − VG FSM: baseline has, VG dropped) | **109** |
| **VG wins** (VG FSM − baseline FSM: VG recovers, baseline misses) | **108** |
| **net VG − baseline** | **−1** |

**VG mode is net-NEUTRAL vs baseline.** Its 108 real copy/multimapper recoveries are almost
exactly cancelled by 109 regressions. The `VG⊇baseline` invariant (UnionBaseline rescue) is
NOT holding — VG leaks 109 baseline-recoverable annotated isoforms.

## Within the st_only gap

Of the 668 st_only misses (ST has, VG missed), **94 (14%) are recovered by baseline `-L`** =
VG dropped a baseline-recoverable, ST-confirmed isoform. Concentrated in no_overlap (47),
partial_multi (27), altsplice (11). Confirmed concretely: NC_073224.2 XM_063708549.1 (46 kb) —
VG emits NOTHING at the locus; baseline emits 8 transcripts incl. the FSM match.

## Why this is the lever

- **Pure recall recovery, zero precision cost:** baseline already produces these 109 correctly;
  VG just needs to stop losing them. No indistinguishability tradeoff (unlike every flow/filter
  recall lever).
- **Unlocks VG's actual value:** fixing the regressions takes VG from −1 to **+108 net** over
  baseline — the copy-finding gain VG is supposed to deliver, currently masked by the leak.
- **Directly answers "extend to VG for more isoforms."**

## Caveat / next step (scoping)

The VG run carried extra precision flags (RUSTLE_VG_DECISIVE_GATE=1, TANDEM). Some of the 109
regressions may be the decisive gate / family-graph / apportionment over-suppressing real
baseline-recoverable transcripts (the gate trades recall for phantom-suppression). The next
step is to attribute the 109 to their VG-pipeline cause:
  - re-run VG with vs without each suppression flag on the regression loci, OR
  - check whether the UnionBaseline rescue (meant to guarantee ⊇baseline) fires for them.
Then fix the dominant cause. Tooling: this analysis reuses `/tmp/strand_safety/all/*/off_fsm.txt`
(baseline) and `perchrom/*/r_fsm.txt` (VG); regression set is `baseline_off − VG` per chrom.
