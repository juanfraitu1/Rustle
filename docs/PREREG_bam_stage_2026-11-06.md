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

## Result (2026-11-06 10:55): panel fetch 189 s, identical; premise corrected; genome-wide arm follows

hot12 (job 4447696): copy 292 s; chr13 fetch from NVMe 189.3 s (bar <= 400 PASS), wall 1,655 s, five tables
identical to h9out/auto. Correction to the premise: the 779 s panel fetch was hot11's (the unit arm), while
hot10's unset arm fetched the same region from NFS in 188.7 s. So NVMe equals the NFS best case on this panel;
what NFS adds is variance (189-868 s on one region across five runs), consistent with contention on the
shared home filesystem rather than a slow path. n = 1 for NVMe; the panel cannot show the variance is gone.

## Amendment 1 (2026-11-06 11:00, before gw12): genome-wide arm with the staged BAM
gw12: singlenode, gw7 recipe exactly (region threads 8, family 4 x star 16, budget 4M, --threads 64), bin10
(flag unset, verified identical on the panel), BAM from the NVMe stage. Bars: assign <= 5,400 s, six tables
identical to gw_off, RSS < 250 GB; fetch sum (936 regions) reported against gw7's 11,308 s and e2e's 8,152 s.
Prediction: fetch sum 5,000-8,000 s; chr13 obtains its budget earlier because the other fetches finish
sooner; wall 5,200-6,000 s, so the bar is a coin flip; identical. Wall limit 2.5 h.
