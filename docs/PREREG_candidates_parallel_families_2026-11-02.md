# PREREG CANDIDATES PARALLEL FAMILIES (register 1506, 2026-11-02)

Why: `candidate_copies` phase 1 (clusters, consensus, refinement, merge; `src/candidate_copies.rs` `for (fi, fam) in families`)
runs the families one after another; each family's work is independent (own tmp dir `f<fi>`, deterministic `sample_net`, the
minimap2 result cache commits atomically). On the r1504 human panel (60 families) the stage ran > 13 min of wall.

Change: phase 1 over a rayon pool of `--threads` / per-family minimap2 threads workers; per-family logs collected and printed in
family order after the pool; `Mm2.calls` becomes atomic. Phases 2-4 unchanged.

Decision rule (pre-declared): ships as default only if, on the r1504 human panel (`/anvme/workspace/v138da13-rustle/nont2t/human`,
same binary flags):
1. `x.cand.candidates.tsv`, `x.cand.contigs.fa`, `x.cand.counts.tsv` and the cluster members table are byte-identical to the r1504
   run; the log's per-family lines are the same multiset;
2. the stage's `clusters` time falls by >= 3x.
If 1 fails: fix or drop (identical products are the contract). If only 2 fails: keep (costs nothing), record the measured gain.
