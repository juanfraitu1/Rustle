# Contiguous-core gate: downstream effect AT SCALE on the biggest-drop contig

Date: 2026-06-16. Confirms the (robustness-fixed) contiguous-core family-merge
gate's REAL downstream output effect on **NC_073228.2** — the chromosome with
the most would-be drops (54 of 135 cross-copy merge pairs have core_cov < 0.13,
40% drop rate; the biggest in the genome-wide measurement `core_gate_gw.md`).
Gate OFF vs ON inside rustle's actual `--vg --vg-layer2` pipeline.

This is the decisive companion to `core_gate_pipeline.md` (which ran the same
test on NC_073235.2 = only 5/67 drops). Question: does 11x more dropped merges
translate into 11x more downstream churn, or worse, lost real transcripts?

## Setup (OOM-safe, reproducible)

- Binary: `/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle`, built
  2026-06-16 16:30 from HEAD `9452964` (has the gate + the robust-POA fix).
- Contig: **NC_073228.2** (gorilla autosome, 195.3 Mb, **371,526 reads** — the
  largest contig; more reads than NC_073235.2's 188k).
- Slices from `/home/juanfra/winloci_scratch/GGO.{bam,fasta}`; guide
  `/tmp/gw/st_NC_073228.2.gtf` (StringTie -L, 4777 tx). `RAYON_NUM_THREADS=4`,
  serial, single process (~0.7 GB peak RSS observed; well under the 17 GB free).
- Gate env: `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE` (default 0.0 = OFF). ON value: **0.13**
  (the bimodal-valley threshold from `core_gate_gw.md`).
- `RUSTLE_VG_CORE_GATE_TRACE=1` on both runs (additive, never changes decisions).

Commands:
```
samtools faidx GGO.fasta NC_073228.2 > NC_073228.2.fasta; samtools faidx NC_073228.2.fasta
samtools view -b GGO.bam NC_073228.2 > NC_073228.2.bam; samtools index NC_073228.2.bam
# OFF
RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 rustle --vg --vg-layer2 \
  --genome-fasta NC_073228.2.fasta -G st_NC_073228.2.gtf -L NC_073228.2.bam -o cgon_off.gtf
# ON @0.13
RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13 \
  rustle --vg --vg-layer2 --genome-fasta NC_073228.2.fasta -G st_NC_073228.2.gtf \
  -L NC_073228.2.bam -o cgon_on.gtf
```

## TL;DR

The gate **fires on exactly 54 merge pairs** (the predicted count), yet the
downstream output change is **smaller than on NC_073235.2** despite 11x more
dropped merges. **Zero transcripts gained or lost** (4989 = 4989). The ONLY real
change is at **one locus (RSTL.262)**: two transcripts lose their spurious
`copy_status "novel"` + `rescue_class "strand_pure_minority"` tags while keeping
identical coordinates and identical abundance. **Recall and precision vs the
reference are byte-identical** (FSM 1733 = 1733; transcript Sn/Pr 23.3/34.7 both
ways). This is benign cleaning — the gate's intended precision behaviour. **No
real isoform is lost.**

## Result 1 — the gate is active: 54 merges dropped

- Fresh OFF trace: 135 Jaccard-passing cross-copy singleton pairs;
  **54 have core_cov < 0.13** (52 deep-lobe < 0.05, 2 in the [0.05, 0.13)
  valley) — reproduces `core_gate_gw.md` exactly.
- ON @0.13: `would_gate=true` count = **54** — fires on precisely those pairs.
  The 81 true-copy pairs (core_cov >= 0.13; 48 in [0.31, 0.95), 29 >= 0.95) are
  untouched. Direction exactly as designed: drop short-block domain-sharers,
  keep real copies.

## Result 2 — downstream output delta (full contig)

| metric | OFF | ON | delta |
|---|---|---|---|
| transcripts | 4989 | 4989 | **0** |
| genes | 1497 | 1497 | 0 |
| `family_id`-tagged tx | 49 | 49 | 0 |
| `copy_id`-tagged tx | 49 | 49 | 0 |
| `copy_status "novel"` | 111 | 109 | **-2** |
| `rescue_class "strand_pure_minority"` | 62 | 60 | **-2** |
| FSM (=) vs reference | 1733 | 1733 | **0** |
| transcript Sn / Pr | 23.3 / 34.7 | 23.3 / 34.7 | 0 / 0 |

Method note: the raw line diff is large (~57k lines) but that is GTF
line-ORDER nondeterminism (rayon), not content. The numbers above come from a
**coordinate-keyed diff** (`/tmp/gw/diff_keyed.py`): each transcript keyed by
(strand, exon-chain). Keyed result:
- chains only in OFF (lost): **0**
- chains only in ON (gained): **0**
- chains in both: **4989**
- chains with a REAL attr change (copy_status / rescue_class / copy_id /
  source / reference_id / family_size): **2**
- chains with ONLY a cosmetic family_id renumber: **0**

So the 4989 = 4989 is NOT cancellation — it is literally the same transcript set
by coordinate, with 2 attribute edits.

## Result 3 — the one real change: RSTL.262 (48,724,206-48,742,369, + strand)

Same 5 transcripts in both runs, **identical exon coordinates and identical
cov / FPKM / TPM / abundance_min** values. The only difference: the two
transcripts sourced from `guide:STRG.440.1` and `guide:STRG.440.5` carry
`rescue_class "strand_pure_minority"` + `copy_status "novel"` with the gate OFF;
with the gate ON those two tags are **dropped**. (Transcript_id .1-.5 ordering
reshuffles, but the set keyed by source/coords/cov is identical.)

Interpretation: removing the spurious cross-copy domain-sharer merges that were
contaminating this family's variation graph cleaned up its false novel-copy
attribution. The transcripts are real and guide-backed (STRG.440.x) — they
**stay emitted**; only the unwarranted "novel copy" / rescue tags go away. This
is the exact analogue of the NC_073235.2 RSTL.647 result in `core_gate_pipeline.md`.

## Verdict

- **Benign-cleaning (expected), not harmful.** 0 transcripts lost, 0 real
  isoforms dropped (FSM 1733 = 1733; Sn/Pr unchanged). The only output change is
  removing 2 spurious `copy_status "novel"` / `rescue_class` tags at one locus —
  a precision improvement on copy attribution, not a recall cost.
- **Drop count does NOT scale to output churn.** 54 dropped merges (11x the 5 on
  NC_073235.2) produce FEWER downstream edits (2 tag cleanings, 0 family
  renumbers) than NC_073235.2 (1 locus cleaned + 2 cosmetic renumbers). The 54
  gated merges sit on family graphs whose final flow-extracted transcripts did
  not depend on those specific cross-copy exon fusions — confirming the gate
  prunes contaminating edges without disturbing the emitted assembly.
- **Honest caveat — the at-scale downstream effect is essentially negligible in
  magnitude** (2 tags on 1 of 1497 genes), even on the worst-case contig. The
  value of the gate is precision hygiene (it removes a real false novel-copy
  call), not a headline metric move. It is also strictly safe: nothing real is
  lost. Per-chrom / within-chromosome families only (OOM protocol); cross-chrom
  paralogs are not the gate's target.

## Artifacts

- GTFs: `/tmp/gw/cgon_off.gtf`, `/tmp/gw/cgon_on.gtf`
- Logs (with CORE_TRACE): `/tmp/gw/cgon_off.log`, `/tmp/gw/cgon_on.log`
- Keyed diff script: `/tmp/gw/diff_keyed.py`
- gffcompare: `/tmp/gw/gc_off.*`, `/tmp/gw/gc_on.*`
- BAM/FASTA slices removed after the run (regenerate via the commands above).
