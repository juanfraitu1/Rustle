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

## ATTRIBUTION (2026-06-09) — cause = core --vg family processing, NOT a tunable gate

Two-stage attribution (`vg_regression_attribution.py` slices + whole-chrom ±flag on NC_073247.2):

**Slice stage (per-locus ±flag, n=109):** **102/109 (94%) are `context_only`** — full VG on the
locus SLICE *recovers* the transcript. So the drop is NOT a per-locus gate (those reproduce on a
slice); it's a **whole-chromosome-scale effect**. On a slice the paralog family doesn't form, so
VG assembles the gene normally; on the full chromosome the locus is absorbed into a family.
(core_vg 5, tandem 1, single_exon 1, decisive_gate 0.)

**Whole-chromosome stage (NC_073247.2, 17 regressions):** removing `RUSTLE_VG_DECISIVE_GATE`
recovers only **1 of the 17** regressions. It reaches baseline's FSM *count* (full 1025 →
no-decisive 1036 = baseline 1036, +11) but via a **set reshuffle** — the +11 are mostly different
(VG-specific) transcripts, not the regressions. So **the decisive gate is NOT the cause** (and
relaxing it just trades for phantoms).

**Conclusion:** the 109 regressions are **core `--vg` family processing** — when a gene is pulled
into a paralog family at whole-chromosome scale, the family apportionment/merge/scope drops it.
This is NOT a removable flag; the fix is in VG family construction/apportionment (connects to the
`vg/flow-capacity-apportionment` branch + the O5/borrow work). Harder than a gate flip, but the
cause is identified and localized to the family stage.

Tooling: `vg_regression_attribution.py`, `vg_regressions.json`. Baseline FSM
`/tmp/strand_safety/all/*/off_fsm.txt`, VG FSM `perchrom/*/r_fsm.txt`.
