# PREREG row 1547: BAM on NVMe workspace instead of NFS home (2026-11-06 10:15, before the copy and run)

## Question
`reads_in_region` on chr13 (1.42M reads) took 269-868 s across identical-code runs with the BAM on the atuin
NFS home; whole-run fetch sums 8,152-11,308 s over 936 regions. Does reading the same BAM from the /anvme
NVMe workspace cut and stabilise the fetch? Identity-preserving by construction (same bytes).

## Arms (hot9 panel, singlenode, hot10 recipe, bin10, flag unset, wall limit 1 h)
- nfs: the hot10 unset arm as measured (fetch 779 s, wall 1,654 s) and hot11 (fetch 779 s)
- nvme: same command with `--bam /anvme/workspace/v138da13-rustle/stage/A119b.t2t.bam` (cp of BAM + .bai)
## Bars
- chr13 fetch <= 400 s (about half of the panel's 779 s) and five tables identical to h9out/auto.
- Copy time recorded (one-off per sample; counts against end to end if adopted).
## Prediction
fetch 200-350 s, wall ~1,100-1,250 s, identical. If it holds, the genome-wide arm is gw7's recipe with the
staged BAM; the pipeline change would be a documented staging step, not code.
